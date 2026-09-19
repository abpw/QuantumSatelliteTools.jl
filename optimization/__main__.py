import argparse
from pathlib import Path

from . import Bayesian, Genetic, JuliaObjective, Problem


def main():
    parser = argparse.ArgumentParser()
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--layout", choices=("population", "random", "equispaced", "grid", "census_grid"), default="population")
    common.add_argument("--region", choices=("us", "china", "global"), default="us", help="census grid region")
    common.add_argument("--ground-stations", type=int, help="ignored for grids; defaults to 101 for equispaced, otherwise 100")
    common.add_argument("--orbits", "--shells", dest="shells", type=int,
                        choices=(1, 2, 3, 5), default=5)
    common.add_argument("--grid-distance", type=int, default=400)
    common.add_argument("--duration", type=float, default=24 * 3600)
    common.add_argument("--step", type=float, default=30)
    common.add_argument("--attempt-rate", type=float, default=1e9)
    common.add_argument("--seed", type=int, default=42)
    commands = parser.add_subparsers(dest="algorithm", required=True)
    bayesian = commands.add_parser("bayesian", parents=[common])
    bayesian.add_argument("--calls", type=int, default=1200)
    bayesian.add_argument("--acquisition", choices=("EI", "LCB", "PI"), default="LCB")
    bayesian.add_argument("--checkpoint", type=Path)
    genetic = commands.add_parser("genetic", parents=[common])
    genetic.add_argument("--calls", type=int, default=1200)
    args = parser.parse_args()
    problem = Problem(shells=args.shells)
    ground_stations = (args.ground_stations if args.ground_stations is not None
                       else 101 if args.layout == "equispaced" else 100)
    objective = JuliaObjective(layout=args.layout, ground_stations=ground_stations,
                               grid_distance_km=args.grid_distance, duration_s=args.duration,
                               step_s=args.step, attempt_rate_hz=args.attempt_rate, seed=args.seed,
                               region=args.region)
    optimizer = (Bayesian(calls=args.calls, acquisition=args.acquisition, seed=args.seed,
                          checkpoint=args.checkpoint) if args.algorithm == "bayesian"
                 else Genetic(calls=args.calls, seed=args.seed))
    result = optimizer.optimize(problem, objective)
    print(result)


if __name__ == "__main__":
    main()
