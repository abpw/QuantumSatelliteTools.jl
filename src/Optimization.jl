using SatelliteToolboxPropagators: Propagators
using ..AstronomyGeometry: seconds_per_day
using ..GenerateGroundStations: generate_population_center_gses, generate_random_gses,
    generate_equispaced_gses, generate_grid_gses, generate_census_grid_gses
using ..GenerateSatellites: generate_regular_array_TLEs
using ..LossCalculation: total_loss
using ..Simulator: Node, NodeLabel, Propagator, satellite, ground_station,
    downlink_pairwise_fn, simulate

struct Shell
    sats_per_plane::Int
    orbital_planes::Int
    inclination_deg::Float64
end

abstract type GroundStationLayout end

struct PopulationLayout <: GroundStationLayout
    count::Int
end

struct RandomLayout <: GroundStationLayout
    count::Int
    seed::Int
    force_land::Bool
end

RandomLayout(count::Int, seed::Int=42) = RandomLayout(count, seed, true)

struct EquispacedLayout <: GroundStationLayout
    count::Int
end

struct GridLayout{F} <: GroundStationLayout
    equatorial_distance_km::Int
    discount_func::F
    force_land::Bool
end

GridLayout(distance::Int=400; discount_func=lat_rad->0, force_land::Bool=false) =
    GridLayout(distance, discount_func, force_land)

"""Precompute the census grid once for repeated constellation evaluations."""
struct CensusGridLayout <: GroundStationLayout
    equatorial_distance_km::Float64
    stations::Vector{Tuple{Float64,Float64}}
end

CensusGridLayout(distance::Real=400, region::AbstractString="us") = CensusGridLayout(
    Float64(distance), generate_census_grid_gses(equatorial_distance_km=distance, region=region))

generate_ground_stations(layout::CensusGridLayout) = copy(layout.stations)

generate_ground_stations(layout::PopulationLayout) = generate_population_center_gses(layout.count)
generate_ground_stations(layout::RandomLayout) = generate_random_gses(
    layout.count; seed=layout.seed, force_land=layout.force_land)
generate_ground_stations(layout::EquispacedLayout) = generate_equispaced_gses(layout.count)
generate_ground_stations(layout::GridLayout) = generate_grid_gses(
    equatorial_distance_km=layout.equatorial_distance_km,
    discount_func=lat_rad -> layout.discount_func(lat_rad),
    force_land=layout.force_land,
)
generate_ground_stations(gses::AbstractVector) = gses
generate_ground_stations(generate::Function) = generate()

abstract type Matching end
struct Unconstrained <: Matching end
struct Constrained <: Matching end

struct DualDownlink{M<:Matching}
    matching::M
end

DualDownlink(; constrained::Bool=true) = DualDownlink(constrained ? Constrained() : Unconstrained())

function link_transmissivities(links)
    by_satellite = Dict{Int,Vector{Tuple{Int,Float64}}}()
    for link in links
        sat, gs = link.node1.label.type == satellite ? (link.node1, link.node2) : (link.node2, link.node1)
        push!(get!(by_satellite, sat.id, Tuple{Int,Float64}[]), (gs.id, 1 - total_loss(link.channel)))
    end
    return by_satellite
end

function (evaluation::DualDownlink{Unconstrained})(links)
    return sum((a[2] * b[2] for channels in values(link_transmissivities(links))
                    for (i, a) in enumerate(channels) for b in channels[i+1:end]); init=0.0)
end

function (evaluation::DualDownlink{Constrained})(links)
    candidates = [(a[2] * b[2], sat, a[1], b[1])
                  for (sat, channels) in link_transmissivities(links)
                  for (i, a) in enumerate(channels) for b in channels[i+1:end]]
    used_satellites, used_gses, score = Set{Int}(), Set{Int}(), 0.0
    for (transmissivity, sat, gs1, gs2) in sort!(candidates; rev=true)
        if sat ∉ used_satellites && gs1 ∉ used_gses && gs2 ∉ used_gses
            score += transmissivity
            push!(used_satellites, sat)
            push!(used_gses, gs1, gs2)
        end
    end
    return score
end

function build_constellation(shells, altitude_km)
    propagators = Propagator[]
    for (i, shell) in enumerate(shells)
        shell.sats_per_plane == 0 && continue
        tles = generate_regular_array_TLEs(
            name_prefix="SHELL $i SAT",
            orbital_planes=shell.orbital_planes,
            sats_per_plane=shell.sats_per_plane,
            altitude_km=altitude_km,
            inclination_rad=deg2rad(shell.inclination_deg),
        )
        append!(propagators, Propagators.init.(Val(:SGP4), tles))
    end
    return propagators
end

"""
Average dual-downlink delivery rate using channel transmissivities and the source
attempt rate. The default matcher is global greedy, not GreedyBackoff.
Assumes one emitted pair per attempt, without SPDC or background-noise modeling.
"""
function evaluate_constellation(
    shells,
    layout;
    duration_s::Real=seconds_per_day,
    step_s::Real=30,
    altitude_km::Int=550,
    constrained::Bool=true,
    evaluation=DualDownlink(constrained=constrained),
    tx_aperture_m::Float64=0.6,
    rx_aperture_m::Float64=0.6,
    wavelength_nm::Float64=1550.0,
    attempt_rate_hz::Real=1e9,
)
    duration_s >= 0 || throw(ArgumentError("duration_s must be nonnegative"))
    step_s > 0 || throw(ArgumentError("step_s must be positive"))
    attempt_rate_hz >= 0 || throw(ArgumentError("attempt_rate_hz must be nonnegative"))
    all(shell -> shell.sats_per_plane >= 0 && shell.orbital_planes > 0 &&
        0 <= shell.inclination_deg <= 180, shells) || throw(ArgumentError("invalid shell"))
    satellites = build_constellation(shells, altitude_km)
    isempty(satellites) && return 0.0
    nodes = vcat(
        [Node(sat, 0.0, NodeLabel(satellite, [])) for sat in satellites],
        [Node(gs, 0.0, NodeLabel(ground_station, [])) for gs in generate_ground_stations(layout)],
    )
    last_sample_s = step_s * max(0, ceil(Int, duration_s / step_s) - 1)
    metrics = simulate(
        nodes, last_sample_s, step_s;
        pairwise_fns=[downlink_pairwise_fn],
        evaluation_fn=links -> evaluation(links),
        tx_aperture_m=tx_aperture_m,
        rx_aperture_m=rx_aperture_m,
        wavelength_nm=wavelength_nm,
        verbose=false,
    )
    return sum(metrics) / length(metrics) * attempt_rate_hz
end

function evaluate_constellation(sats_per_plane, orbital_planes, inclinations_deg, layout; kwargs...)
    length(sats_per_plane) == length(orbital_planes) == length(inclinations_deg) ||
        throw(ArgumentError("shell parameters must have equal lengths"))
    shells = Shell.(Int.(sats_per_plane), Int.(orbital_planes), Float64.(inclinations_deg))
    return evaluate_constellation(shells, layout; kwargs...)
end
