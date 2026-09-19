from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np

_RUNTIME = None


@dataclass(frozen=True)
class Problem:
    shells: int = 5
    total_satellites: int = 100
    min_satellites: int = 0

    def __post_init__(self):
        if (self.shells < 1 or self.total_satellites < 0 or self.min_satellites < 0
                or self.min_satellites * self.shells > self.total_satellites):
            raise ValueError("invalid constellation bounds")

    @property
    def bounds(self):
        return [(-5.0, 5.0)] * (self.shells - 1) + [(0.0, 180.0)] * self.shells

    def decode(self, parameters):
        y = np.r_[parameters[: self.shells - 1], 0.0]
        proportions = np.exp(y - y.max())
        proportions /= proportions.sum()
        remaining = self.total_satellites - self.shells * self.min_satellites
        raw = self.min_satellites + proportions * remaining
        satellites = np.floor(raw).astype(int)
        order = np.argsort(-(raw - satellites), kind="stable")
        satellites[order[: self.total_satellites - satellites.sum()]] += 1
        inclinations = parameters[self.shells - 1 :]
        return np.column_stack((satellites, np.ones(self.shells, dtype=int), inclinations))


@dataclass(frozen=True)
class JuliaObjective:
    layout: str = "population"
    ground_stations: int = 100
    grid_distance_km: int = 400
    seed: int = 42
    duration_s: float = 24 * 3600
    step_s: float = 30
    constrained: bool = True
    attempt_rate_hz: float = 1e9
    region: str = "us"

    def __post_init__(self):
        if (self.layout not in {"population", "random", "equispaced", "grid", "census_grid"}
                or self.region not in {"us", "china", "global"}
                or self.ground_stations < 0 or self.grid_distance_km <= 0
                or self.duration_s < 0 or self.step_s <= 0
                or self.attempt_rate_hz < 0
                or self.layout == "equispaced" and self.ground_stations % 2 == 0):
            raise ValueError("invalid Julia objective settings")

    @cached_property
    def _census_grid(self):
        return _julia().CensusGridLayout(self.grid_distance_km, self.region)

    def __call__(self, configuration):
        values = np.asarray(configuration, dtype=float)
        if (values.ndim != 2 or values.shape[1] != 3 or not np.isfinite(values).all()
                or (values[:, :2] != np.floor(values[:, :2])).any()):
            raise ValueError("configuration must contain integer counts and finite inclinations")
        optimization = _julia()
        layouts = {
            "population": lambda: optimization.PopulationLayout(self.ground_stations),
            "random": lambda: optimization.RandomLayout(self.ground_stations, self.seed),
            "equispaced": lambda: optimization.EquispacedLayout(self.ground_stations),
            "grid": lambda: optimization.GridLayout(self.grid_distance_km),
            "census_grid": lambda: self._census_grid,
        }
        layout = layouts[self.layout]()
        return float(optimization.evaluate_constellation(
            values[:, 0].astype(int), values[:, 1].astype(int), values[:, 2], layout,
            duration_s=self.duration_s, step_s=self.step_s, constrained=self.constrained,
            attempt_rate_hz=self.attempt_rate_hz,
        ))


def _julia():
    global _RUNTIME
    if _RUNTIME is None:
        from juliacall import Main as jl

        jl.seval("using QuantumSatelliteTools")
        _RUNTIME = jl.QuantumSatelliteTools.Optimization
    return _RUNTIME


@dataclass(frozen=True)
class Bayesian:
    calls: int = 1200
    initial_points: int = 100
    acquisition: str = "LCB"
    xi: float = 0.0
    seed: int = 42
    checkpoint: Path | None = None

    def __post_init__(self):
        if self.calls < 1 or self.initial_points < 0 or not np.isfinite(self.xi) or self.xi < 0:
            raise ValueError("invalid Bayesian optimizer settings")

    def optimize(self, problem, objective):
        from skopt import gp_minimize, load
        from skopt.callbacks import CheckpointSaver
        from skopt.space import Real

        metadata = self.checkpoint.with_suffix(self.checkpoint.suffix + ".txt") if self.checkpoint else None
        signature = repr((problem, objective, self.acquisition, self.seed, self.initial_points, self.xi))
        if self.checkpoint and self.checkpoint.exists() and (
                not metadata.exists() or metadata.read_text() != signature):
            raise ValueError("checkpoint belongs to another experiment")
        if metadata:
            metadata.parent.mkdir(parents=True, exist_ok=True)
            metadata.write_text(signature)
        previous = load(self.checkpoint) if self.checkpoint and self.checkpoint.exists() else None
        used = len(previous.x_iters) if previous else 0
        if used >= self.calls:
            return previous
        callback = [CheckpointSaver(self.checkpoint, store_objective=False)] if self.checkpoint else []
        return gp_minimize(
            lambda parameters: -objective(problem.decode(parameters)),
            [Real(low, high) for low, high in problem.bounds],
            n_calls=self.calls - used,
            n_initial_points=min(max(0, self.initial_points - used), self.calls - used),
            acq_func=self.acquisition,
            xi=self.xi,
            x0=previous.x_iters if previous else None,
            y0=previous.func_vals if previous else None,
            random_state=previous.random_state if previous else self.seed,
            callback=callback,
        )


@dataclass(frozen=True)
class GeneticResult:
    parameters: np.ndarray
    configuration: np.ndarray
    reward: float


@dataclass(frozen=True)
class Genetic:
    calls: int = 1200
    population: int = 25
    mutation_rate: float = 0.5
    mutation_decay: float = 0.01
    minimum_mutation_rate: float = 0.05
    parents: int = 5
    random_fraction: float = 0.1
    seed: int = 42

    def __post_init__(self):
        if (self.calls < 1 or self.population < self.parents or self.parents < 2
                or self.mutation_decay < 0
                or not 0 <= self.minimum_mutation_rate <= self.mutation_rate <= 1
                or not 0 <= self.random_fraction <= 1):
            raise ValueError("invalid genetic optimizer settings")

    def optimize(self, problem, objective):
        rng = np.random.default_rng(self.seed)
        bounds = np.asarray(problem.bounds)
        population = rng.uniform(bounds[:, 0], bounds[:, 1], (self.population, len(bounds)))
        best = None
        for generation, used in enumerate(range(0, self.calls, self.population)):
            batch = population[: min(self.population, self.calls - used)]
            rewards = np.asarray([objective(problem.decode(parameters)) for parameters in batch])
            winner = int(rewards.argmax())
            if best is None or rewards[winner] > best.reward:
                best = GeneticResult(batch[winner].copy(), problem.decode(batch[winner]), float(rewards[winner]))
            if used + len(batch) == self.calls:
                break
            parents = batch[np.argsort(rewards)[-self.parents :]]
            offspring = int(np.ceil(self.population * (1 - self.random_fraction)))
            population = rng.uniform(bounds[:, 0], bounds[:, 1], (self.population, len(bounds)))
            mutation_rate = max(self.mutation_rate * np.exp(-self.mutation_decay * generation),
                                self.minimum_mutation_rate)
            for i in range(offspring):
                first, second = parents[rng.choice(len(parents), 2, replace=False)]
                point = rng.integers(1, len(bounds)) if len(bounds) > 1 else 1
                child = np.r_[first[:point], second[point:]]
                mutate = rng.random(len(bounds)) < mutation_rate
                child[mutate] = rng.uniform(bounds[mutate, 0], bounds[mutate, 1])
                population[i] = child
        return best


__all__ = ["Bayesian", "Genetic", "GeneticResult", "JuliaObjective", "Problem"]
