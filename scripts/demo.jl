using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf: @printf
using Statistics: mean

using SatelliteToolboxPropagators: Propagators

using QuantumSatelliteTools
using QuantumSatelliteTools.constants: SECONDS_PER_DAY
using QuantumSatelliteTools.types: GS, Propagator
using QuantumSatelliteTools.helpers: _to_float64
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel
using QuantumSatelliteTools.GenerateGroundStations: generate_equispaced_gses
using QuantumSatelliteTools.GenerateSatellites: generate_regular_array_TLEs
using QuantumSatelliteTools.LossCalculation: total_loss

"""
    build_constellation(
        name_prefix::String,
        orbital_planes::Integer,
        sats_per_plane::Integer, 
        altitude_km::Real, 
        inclination_rad::Real, 
    ) -> Vector{Propagator}

Build simple synthetic satellite constellation and initialize SGP4 propagators.

# Arguments

  - `name_prefix::String`: Satellite naming prefix.
  - `orbital_planes::Integer`: Number of orbital planes.
  - `sats_per_plane::Integer`: Satellites per plane.
  - `altitude_km::Real`: Orbital altitude [km].
  - `inclination_rad::Real`: Orbital inclination [rad].

# Returns

  - `Vector{Propagator}`: SGP4 propagators for generated satellites.
"""
function build_constellation(;
    name_prefix::String,
    orbital_planes::Integer,
    sats_per_plane::Integer,
    altitude_km::Real,
    inclination_rad::Real,
)
    tles = generate_regular_array_TLEs(;
        name_prefix=name_prefix,
        orbital_planes=orbital_planes,
        sats_per_plane=sats_per_plane,
        altitude_km=altitude_km,
        inclination_rad=inclination_rad,
    )
    return [Propagators.init(Val(:SGP4), tle) for tle in tles]
end

"""
    evaluate_constellation(
        propagators::Vector{Propagator},
        gses::Vector{GS};
        duration_s::Integer=5400,
        step_s::Integer=300
    ) -> NamedTuple

Evaluate constellation performance metrics on a grid of times and ground stations.

Metrics:

  - `coverage_fraction`: fraction of (station, time) points with >= 1 visible satellite.
  - `mean_visible_satellites`: average number of visible satellites per (station, time).
  - `mean_best_link_transmissivity`: mean of the best visible link transmissivity per
    covered point.

# Arguments

  - `propagators::Vector{Propagator}`: Satellite propagators.
  - `gses::Vector{GS}`: Ground station locations.
  - `duration_s::Real=5400`: Simulation duration [s].
  - `step_s::Real=300`: Time step [s].

# Returns

  - `NamedTuple` with fields:
      + `satellites::Integer`
      + `coverage_fraction::Real`
      + `mean_visible_satellites::Real`
      + `mean_best_link_transmissivity::Real`
"""
function evaluate_constellation(
    propagators::Vector{Propagator},
    gses::Vector{GS};
    duration_s::Real=5400,
    step_s::Real=300,
)
    t₀ = propagators[1].sgp4d.epoch

    total_points = 0
    covered_points = 0
    visible_sat_sum = 0
    best_link_transmissivities = Float64[]

    for Δt_s in 0:step_s:duration_s
        t = t₀ + Δt_s / SECONDS_PER_DAY

        for gs in gses
            total_points += 1
            visible_count = 0
            best_transmissivity = 0.0

            for sat in propagators
                channel = FreespaceChannel(sat, gs; time=t)
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

    mean_best_link_transmissivity =
        isempty(best_link_transmissivities) ? 0.0 : mean(best_link_transmissivities)

    return (
        satellites=length(propagators),
        coverage_fraction=covered_points / total_points,
        mean_visible_satellites=visible_sat_sum / total_points,
        mean_best_link_transmissivity=mean_best_link_transmissivity,
    )
end

"""
    print_comparison(
        title_a::String, metrics_a::NamedTuple, title_b::String, metrics_b::NamedTuple
    ) -> Nothing

Print a compact comparison table of two constellation configurations.

# Arguments

  - `title_a::String`: First configuration title.
  - `metrics_a::NamedTuple`: First configuration metrics.
  - `title_b::String`: Second configuration title.
  - `metrics_b::NamedTuple`: Second configuration metrics.
"""
function print_comparison(
    title_a::String, metrics_a::NamedTuple, title_b::String, metrics_b::NamedTuple
)
    println("\nConstellation Comparison")
    println("---------------------------------------------------------------")
    @printf("%-32s %-16s %-16s\n", "Metric", title_a, title_b)
    println("---------------------------------------------------------------")
    @printf("%-32s %-16d %-16d\n", "Satellites", metrics_a.satellites, metrics_b.satellites)
    @printf(
        "%-32s %-16.2f %-16.2f\n",
        "Coverage (% of station-time)",
        100 * metrics_a.coverage_fraction,
        100 * metrics_b.coverage_fraction
    )
    @printf(
        "%-32s %-16.2f %-16.2f\n",
        "Mean visible satellites",
        metrics_a.mean_visible_satellites,
        metrics_b.mean_visible_satellites
    )
    @printf(
        "%-32s %-16.4f %-16.4f\n",
        "Mean best-link transmissivity",
        metrics_a.mean_best_link_transmissivity,
        metrics_b.mean_best_link_transmissivity
    )
    println("---------------------------------------------------------------")
end

# 1) Build a small deterministic set of ground stations.
# We use an odd n so the Fibonacci-based equispaced generator is valid.
gses::Vector{GS} = [_to_float64(gs) for gs in generate_equispaced_gses(9)]

# 2) Define two simple synthetic constellations to compare.
# A = sparse, B = denser constellation with more orbital planes/satellites.
constellation_a = build_constellation(;
    name_prefix="SPARSE",
    orbital_planes=4,
    sats_per_plane=4,
    altitude_km=550,
    inclination_rad=deg2rad(53.0),
)
constellation_b = build_constellation(;
    name_prefix="DENSE",
    orbital_planes=8,
    sats_per_plane=8,
    altitude_km=550,
    inclination_rad=deg2rad(53.0),
)

# 3) Simulate both over a short time horizon and compare quantitative metrics.
metrics_a = evaluate_constellation(constellation_a, gses; duration_s=5400, step_s=300)
metrics_b = evaluate_constellation(constellation_b, gses; duration_s=5400, step_s=300)

print_comparison("Sparse (4x4)", metrics_a, "Dense (8x8)", metrics_b)

println("\nInterpretation:")
println(
    "- Higher coverage means more ground stations have " *
        "at least one visible satellite at sampled times.",
)
println(
    "- Higher mean best-link transmissivity means " *
        "the strongest available link is better on average.",
)

