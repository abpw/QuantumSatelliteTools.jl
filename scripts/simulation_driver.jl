using SatelliteToolboxTle
using GenerateGroundStations
using GenerateTLEs
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using AstronomyGeometry: FreespaceChannel, transmissivity, decibels_to_probability
using SimpleWeightedGraphs, Graphs
using LossCalculation: reflector_loss, swapping_loss

GSES = [
    # North America
    (40.7128, -74.0060),   # New York City, USA (1)
    (34.0522, -118.2437),  # Los Angeles, USA (2)
    (19.4326, -99.1332),  # Mexico City, Mexico (3)
    (49.2827, -123.1207),  # Vancouver, Canada
    # South America
    (-23.5505, -46.6333),  # São Paulo, Brazil
    (-34.6037, -58.3816),  # Buenos Aires, Argentina
    (4.7110, -74.0721),  # Bogotá, Colombia
    # Europe
    (51.5074, -0.1278),   # London, UK
    (48.8566, 2.3522),   # Paris, France
    (52.5200, 13.4050),   # Berlin, Germany
    (41.0082, 28.9784),   # Istanbul, Türkiye
    # Africa
    (6.5244, 3.3792),   # Lagos, Nigeria
    (30.0444, 31.2357),   # Cairo, Egypt
    (-26.2041, 28.0473),   # Johannesburg, South Africa
    # Asia
    (39.9042, 116.4074),   # Beijing, China
    (35.6895, 139.6917),   # Tokyo, Japan
    (19.0760, 72.8777),   # Mumbai, India
    (-6.2088, 106.8456),   # Jakarta, Indonesia
    # Oceania
    (-33.8688, 151.2093),  # Sydney, Australia
    (-36.8485, 174.7633)   # Auckland, New Zealand
]

@enum Entities begin
    sat
    gs
    aux_gs
end

"""
Conduct a simulation of a quantum satellite network.

# Arguments
- `duration_s::Int64`: The total duration of the simulation in seconds.
- `interval_s::Int64`: The interval at which the system is sampled in seconds.


# Returns
- An aggregate representation of the performance of the system.
"""
function simulation_driver(duration_s::Int64, interval_s::Int64,
        epoch::Float64; gses::Vector{Tuple{Float64}}=GSES,
        aux_gses::Union{Vector{Tuple{Float64}}, Missing}=missing)
    # Downloading Starlink TLEs 8.21.2025
    tles = read_tles_from_file(joinpath(@__DIR__, "../databases/starlink.tle"))
    propagators = [Propagators.init(Val(:SGP4), tle) for tle in tles]
    aux_gses = ismissing(aux_gses) ? [] : aux_gses

    # Map entity to an integer index
    node_map = Dict{Union{OrbitPropagatorSgp4,Tuple{Float64}},Int64}(
        node => index for (index, node) in enumerate(vcat(propagators, gses, aux_gses)))

    # Map integer index of entity to entity type
    types = Dict{Int64, Entities}(
        node_map[propagator] => sat for propagator in propagators)
    merge!(types, Dict{Int64, Entities}(
        node_map[gs] => gs for gs in gses))
    merge!(types, Dict{Int64, Entities}(
        node_map[aux_gs] => aux_gs for aux_gs in aux_gses))

    for experiment in [Val(:REF), Val(:DD)]
        for time_elapsed in 0:interval_s:duration_s-1
            sources, destinations, weights =
                Vector{Int64}(), Vector{Int64}(), Vector{Float64}()
            for propagator in propagators
                for gs in gses
                    channel = FreespaceChannel(propagator, gs, time=epoch+time_elapsed)
                    if !isnothing(channel)
                        push!(sources, node_map[propagator])
                        push!(destinations, node_map[gs])
                        push!(weights, transmissivity(channel))
                    end
                end
            end
            graph = SimpleWeightedGraph(sources, destinations, weights)
        end
    end
end

function all_pairs_path_probs(g::SimpleWeightedGraph, types::Dict{Int64, Entities}, experiment)
    path_probs = zeros(nv(g), nv(g)) + I
    for edge ∈ edges(g)
        source = min(src(edge), dst(edge))
        destination = max(src(edge), dst(edge))
        path_probs[source, destination] = g.weights[source, destination]
    end
    for k ∈ 2:nv(g)-1
        if experiment == Val(:REF)
            if types[k] != sat
                continue
            end
            node_cost = decibel_to_probability(reflector_loss())
        end
        if experiment == Val(:DD)
            node_cost = types[k] == sat ? 1 : swapping_loss()
        end
        for i ∈ 1:k-1
            if experiment == Val(:DD) && !xor(types[k] == sat, types[i] == sat)
                continue
            end
            for j ∈ k+1:nv(g)
                if experiment == Val(:DD) && !xor(types[k] == sat, types[j] == sat)
                    continue
                end
                path_probs[i, j] = min(path_probs[i, j], path_probs[i, k] * path_probs[k, j] * node_cost)
            end
        end
    end
    return sum([path_probs[i,j] for i ∈ vertices(g) if types[i] == gs for j ∈ vertices(g) if types[j] == gs && j != i])
end