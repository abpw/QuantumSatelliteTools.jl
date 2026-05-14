using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf
using Statistics
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using QuantumSatelliteTools
using QuantumSatelliteTools.GenerateGroundStations: generate_equispaced_gses
using QuantumSatelliteTools.GenerateSatellites: generate_regular_array_TLEs
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel
using QuantumSatelliteTools.LossCalculation: total_loss

# Ground stations are represented as (latitude_rad, longitude_rad) tuples.
const GroundStation = Tuple{Float64, Float64}
const Propagator = OrbitPropagatorSgp4{Float64, Float64}

"""
Build a simple synthetic constellation and initialize SGP4 propagators.
"""
function build_constellation(;
    name_prefix::String,
    orbits::Int,
    sats_per_orbit::Int,
    altitude_km::Int,
    inclination_rad::Float64,
)
    tles = generate_regular_array_TLEs(
        name_prefix=name_prefix,
        orbits=orbits,
        sats_per_orbit=sats_per_orbit,
        altitude_km=altitude_km,
        inclination_rad=inclination_rad,
    )
    return [Propagators.init(Val(:SGP4), tle) for tle in tles]
end

"""
Evaluate constellation performance on a grid of times and ground stations.

Metrics:
- coverage_fraction: fraction of (station, time) points with >= 1 visible satellite.
- mean_visible_satellites: average number of visible satellites per (station, time).
- mean_best_link_transmissivity: mean of the best visible link transmissivity per covered point.
"""
function evaluate_constellation(
    propagators::Vector{Propagator},
    gses::Vector{GroundStation};
    duration_s::Int = 5400,
    step_s::Int = 300,
)
    t0 = propagators[1].sgp4d.epoch

    total_points = 0
    covered_points = 0
    visible_sat_sum = 0
    best_link_transmissivities = Float64[]

    for dt_s in 0:step_s:duration_s
        t = t0 + dt_s / 86400
        for gs in gses
            total_points += 1
            visible_count = 0
            best_transmissivity = 0.0

            for sat in propagators
                channel = FreespaceChannel(sat, gs, time=t)
                if !isnothing(channel)
                    visible_count += 1
                    # total_loss returns failure probability, so success = 1 - loss.
                    link_transmissivity = 1 - total_loss(channel)
                    best_transmissivity = max(best_transmissivity, link_transmissivity)
                end
            end

            visible_sat_sum += visible_count
            if visible_count > 0
                covered_points += 1
                push!(best_link_transmissivities, best_transmissivity)
            end
        end
    end

    return (
        satellites=length(propagators),
        coverage_fraction=covered_points / total_points,
        mean_visible_satellites=visible_sat_sum / total_points,
        mean_best_link_transmissivity=isempty(best_link_transmissivities) ? 0.0 : mean(best_link_transmissivities),
    )
end

"""
Print a compact comparison table for two constellation configurations.
"""
function print_comparison(title_a::String, metrics_a, title_b::String, metrics_b)
    println("\nConstellation Comparison")
    println("---------------------------------------------------------------")
    @printf("%-32s %-16s %-16s\n", "Metric", title_a, title_b)
    println("---------------------------------------------------------------")
    @printf("%-32s %-16d %-16d\n", "Satellites", metrics_a.satellites, metrics_b.satellites)
    @printf("%-32s %-16.2f %-16.2f\n", "Coverage (% of station-time)", 100 * metrics_a.coverage_fraction, 100 * metrics_b.coverage_fraction)
    @printf("%-32s %-16.2f %-16.2f\n", "Mean visible satellites", metrics_a.mean_visible_satellites, metrics_b.mean_visible_satellites)
    @printf("%-32s %-16.4f %-16.4f\n", "Mean best-link transmissivity", metrics_a.mean_best_link_transmissivity, metrics_b.mean_best_link_transmissivity)
    println("---------------------------------------------------------------")
end

# 1) Build a small deterministic set of ground stations.
# We use an odd n so the Fibonacci-based equispaced generator is valid.
gses = GroundStation[(Float64(gs[1]), Float64(gs[2])) for gs in generate_equispaced_gses(9)]

# 2) Define two simple synthetic constellations to compare.
# A = sparse, B = denser constellation with more orbital planes/satellites.
constellation_a = build_constellation(
    name_prefix="SPARSE",
    orbits=4,
    sats_per_orbit=4,
    altitude_km=550,
    inclination_rad=deg2rad(53.0),
)
constellation_b = build_constellation(
    name_prefix="DENSE",
    orbits=8,
    sats_per_orbit=8,
    altitude_km=550,
    inclination_rad=deg2rad(53.0),
)

# 3) Simulate both over a short time horizon and compare quantitative metrics.
metrics_a = evaluate_constellation(constellation_a, gses, duration_s=5400, step_s=300)
metrics_b = evaluate_constellation(constellation_b, gses, duration_s=5400, step_s=300)

print_comparison("Sparse (4x4)", metrics_a, "Dense (8x8)", metrics_b)

println("\nInterpretation:")
println("- Higher coverage means more ground stations have at least one visible satellite at sampled times.")
println("- Higher mean best-link transmissivity means the strongest available link is better on average.")

