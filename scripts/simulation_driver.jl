using SatelliteToolboxTle
using QuantumSatelliteTools.GenerateGroundStations
using QuantumSatelliteTools.GenerateTLEs
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using QuantumSatelliteTools.AstronomyGeometry: GS
using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight
using Graphs: neighbors, vertices, nv
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel
using QuantumSatelliteTools.LossCalculation: reflector_loss, swapping_loss, total_loss
using QuantumSatelliteTools.GenerateGroundStations: generate_population_center_gses
using QuantumSatelliteTools.GenerateTLEs: generate_regular_array_TLEs
using ProgressBars: ProgressBar
using DataStructures: PriorityQueue, enqueue!, dequeue!

@enum Entities begin
    satellite
    ground_station
    aux_ground_station
end

@enum Experiment begin
    reflector
    dual_downlink
end

struct Path
    nodes::Vector{Int64}
    transmissivity::Float64
end

struct PathProbs
    num_nodes::Int64
    transmissivity::Float64
end

mutable struct HopProb
    count::Int64
    sum_transmissivity::Float64
end

Entity = Union{OrbitPropagatorSgp4{Float64, Float64},GS}
Propagator = OrbitPropagatorSgp4{Float64, Float64}
GroundStation = Tuple{Number, Number}

"""
Helper to add edge data to three edge data arrays.

# Arguments
- `source::Union{OrbitPropagatorSgp4{Float64, Float64},GS}`: A source entity.
- `destination::Union{OrbitPropagatorSgp4{Float64, Float64},GS}`: A destination
        entity.
- `node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64}`: A
        node map of entities to integer values.
- `time::Float64`: The time of the transmission.
- `sources::Vector{Int64}`: An array of sources for a graph.
- `destinations::Vector{Int64}`: An array of destinations for a graph.
- `weights::Vector{Float64}`: An array of weights for a graph.

# Returns
- Nothing
"""
function edge_data!(
        source::Entity, destination::Entity,
        node_map::Dict{Entity,Int64}, time::Float64,
        sources::Vector{Int64}, destinations::Vector{Int64},
        weights::Vector{Float64}
)
    channel = FreespaceChannel(source, destination, time=time)

    if !isnothing(channel)
        push!(sources, node_map[source])
        push!(destinations, node_map[destination])
        push!(weights, 1-total_loss(channel))
    end
end

"""
Create a graph representing a quatum satellite network at a point in time.

# Arguments
- `propagators::Vector{OrbitPropagatorSgp4}`: An array of orbit propagators
        representing the satellites in the network.
- `gses::Vector{GS}`: An array of tuples holding the coordinates of the ground
        stations in the network.
- `node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64}`: A map
        of entities to integer values.
- `time::Float64`: The time that the graph takes place.
- `experiment::Experiment`: The type of experiment.

# Returns
- SimpleWeightedGraph: A graph representing a quantum satellite network with
        nodes of type Union{OrbitPropagatorSgp4{Float64, Float64},GS}
        representing satellites and ground stations and edges with weights
        representing channel transmissivity.
"""
function make_graph(
        propagators::Vector{Propagator},
        gses::Vector{GroundStation},
        node_map::Dict{Entity,Int64},
        time::Float64,
        experiment::Experiment
)
    sources = Vector{Int64}()
    destinations = Vector{Int64}()
    weights = Vector{Float64}()

    for propagator ∈ propagators
        for gs ∈ gses
            edge_data!(propagator, gs, node_map, time, sources, destinations,
                       weights)
        end
    end

    if experiment == reflector
        for i ∈ eachindex(propagators)
            for j ∈ i+1:length(propagators)
                edge_data!(propagators[i], propagators[j], node_map, time,
                           sources, destinations, weights)
            end
        end
    end

    return SimpleWeightedGraph(sources, destinations, weights)
end

"""
Write the output of a simulation to three designated files.

# Arguments
- `time_transmissivities::Vector{Tuple{Float64, Float64}}`: An array of pairs of
        time stamps and aggregate transmissivities at the given time.
- `hop_probs::Dict{Int64, HopProb}`
- `experiment_prob::Float64`
- `time_file::Union{String, Missing}`: (default missing)
- `hop_file::Union{String, Missing}`: (default missing)
- `agg_file::Union{String, Missing}`: (default missing)

# Returns
- Nothing 
"""

function write_output(
        time_transmissivities::Vector{Tuple{Float64, Float64}},
        hop_probs::Dict{Int64, Float64},
        experiment_prob::Float64,
        time_file::Union{String, Missing}=missing,
        hop_file::Union{String, Missing}=missing,
        agg_file::Union{String, Missing}=missing
)
    if !ismissing(time_file)
        open(pwd()*"/data/"*time_file, "w") do file
            for transmissivity ∈ time_transmissivities
                write(file, "$(transmissivity[1]),$(transmissivity[2])\n")
            end
        end
    end

    if !ismissing(hop_file)
        open(pwd()*"/data/"*hop_file, "w") do file
            for (num_hops, hop_prob) ∈ hop_probs
                write(file, "$(num_hops),$(hop_prob)\n")
            end
        end
    end

    if !ismissing(agg_file)
        open(pwd()*"/data/"*agg_file, "w") do file
            write(file, "$(experiment_prob)")
        end
    end
end

"""
Generate synthetic satellite TLEs and, if needed, population center ground
stations for simulation.

# Arguments
- ``

# Returns
- 
"""
function generate_entities(
        num_sats:: Int64,
        default_gses::Union{Vector{GroundStation}, Missing}=missing,
        default_aux_gses::Union{Vector{GroundStation}, Missing}=missing
)
    println("Generating TLEs")
    tles::Vector{TLE} = [tle
                         for θ ∈ 0.0:18.0:162.0
                         for tle ∈ generate_regular_array_TLEs(
                            orbits=10,
                            sats_per_orbit=num_sats,
                            inclination_rad=θ)]
    propagators = [Propagators.init(Val(:SGP4), tle) for tle ∈ tles]

    gses::Vector{GroundStation} = (ismissing(default_gses)
                                   ? generate_population_center_gses(100)
                                   : default_gses)

    aux_gses::Vector{GroundStation} = (ismissing(default_aux_gses)
                                       ? []
                                       : default_aux_gses)
    
    return propagators, gses, aux_gses
end

"""
Create maps of entities to integer indices.

# Arguments
- ``

# Returns
- 
"""
function make_mappings(
        propagators::Vector{Propagator},
        gses::Vector{GroundStation},
        aux_gses::Vector{GroundStation}
)
    # Map entity to an integer index
    println("Mapping entities to integers")
    node_map = Dict{Entity, Int64}(
        node => index
        for (index, node) ∈ enumerate(vcat(propagators, gses, aux_gses)))

    # Map integer index of entity to entity type
    println("Mapping node integer values to types")
    types = Dict{Int64, Entities}(
        node_map[propagator] => satellite for propagator ∈ propagators)
    merge!(types, Dict{Int64, Entities}(
        node_map[gs] => ground_station for gs ∈ gses))
    merge!(types, Dict{Int64, Entities}(
        node_map[aux_gs] => aux_ground_station for aux_gs ∈ aux_gses))

    return node_map, types
end

"""
Conduct a simulation of a quantum satellite network.

# Example use
- `simulate(86400, 300, 2460900.0, reflector, 5, time_file="time.csv",
        hop_file="hops.csv", agg_file="agg.csv")`

# Arguments
- `duration_s::Int64`: The total duration of the simulation in seconds.
- `interval_s::Int64`: The interval at which the system is sampled in seconds.
- `epoch::Float64`: The time at which the simulation begins.
- `gses::Vector{GS}`: An array of ground stations (default GSES).
- `aux_gses::Vector{GS}`: An array of auxiliary groundstations as
        produced by the DD experiment (default missing).
- `experiment::Experiment`: Which experiment to run.
- `num_sats::Int64`: The number of satellites per orbit.

# Returns
- An aggregate representation of the performance of the system under the
  architecture dictated by the experiment.
"""
function simulate(
        duration_s::Int64,
        interval_s::Int64,
        epoch::Float64,
        experiment::Experiment,
        num_sats::Int64;
        default_gses::Union{Vector{GroundStation}, Missing}=missing,
        default_aux_gses::Union{Vector{GroundStation}, Missing}=missing,
        time_file::Union{String, Missing}=missing,
        hop_file::Union{String, Missing}=missing,
        agg_file::Union{String, Missing}=missing
)
    propagators, gses, aux_gses = generate_entities(num_sats, default_gses,
                                                    default_aux_gses)
    node_map, types = make_mappings(propagators, gses, aux_gses)

    time_transmissivities::Vector{Tuple{Float64, Float64}} = []

    parent_matrix::Matrix{PathProbs} = [
        PathProbs(0, 0) 
        for _ ∈ 1:length(gses), _ ∈ 1:length(gses)]

    println("Beginning simulation")
    for time_elapsed ∈ ProgressBar(0:interval_s:duration_s-1)
        graph = make_graph(propagators, vcat(gses, aux_gses), node_map,
                           epoch+time_elapsed, experiment)
        all_pairs_path_probs!(graph, types, experiment, parent_matrix)
        push!(time_transmissivities,
             (time_elapsed, aggregate_paths_probs(parent_matrix, length(gses))))
    end

    hop_probs = aggregate_hop_probs(parent_matrix)
    experiment_prob = sum(values(hop_probs)) / length(hop_probs)
    write_output(time_transmissivities, hop_probs, experiment_prob, time_file,
                 hop_file, agg_file)
    
    return experiment_prob
end

function all_pairs_path_probs!(
        g::SimpleWeightedGraph,
        types::Dict{Int64, Entities},
        experiment::Experiment,
        path_probs::Matrix{PathProbs}
)
    gses = [i for i ∈ vertices(g) if types[i] == ground_station]
    gs_idxes = Dict{Int64, Int64}(int => idx for (idx, int) ∈ enumerate(gses))

    for i ∈ 1:length(gses), j ∈ 1:length(gses)
        path_probs[i, j] = PathProbs(0, 0)
    end

    node_cost_sat = experiment == reflector ? 1 - reflector_loss() : 1
    node_cost_gs = experiment == reflector ? 1 : 1 - swapping_loss()

    for gs ∈ gses
        dists::Vector{Float64} = [(i == gs ? 1.0 : 0.0) for i ∈ 1:nv(g)]
        prevs::Vector{Union{UndefInitializer, Int64}} = [undef for _ ∈ 1:nv(g)]
        priority_queue = PriorityQueue()
        seen = Set{Int64}()
        push!(seen, gs)

        for vertex ∈ vertices(g)
            enqueue!(priority_queue, vertex => vertex == gs ? 0 : Inf)
        end

        while length(priority_queue) > 0
            curr = dequeue!(priority_queue)
            push!(seen, curr)
            node_cost = types[curr] == satellite ? node_cost_sat : node_cost_gs

            for neighbor ∈ neighbors(g, curr)
                if !(neighbor ∈ seen)
                    alt_prob = (dists[curr]
                                * node_cost
                                * get_weight(g, curr, neighbor))

                    if alt_prob > dists[neighbor]
                        prevs[neighbor] = curr
                        dists[neighbor] = alt_prob
                        priority_queue[neighbor] = -alt_prob
                    end
                end
            end
        end

        path_probs[gs_idxes[gs], gs_idxes[gs]] = PathProbs(0, 0)

        for gs_end ∈ gses
            if gs_end == gs
                continue
            end

            path_length = 1
            current = gs_end
            continue_flag = false

            while current != gs
                if current == undef
                    path_probs[gs_idxes[gs], gs_idxes[gs_end]] = PathProbs(0, 0)
                    continue_flag = true
                    break
                end

                path_length += 1
                current = prevs[current]
            end

            if continue_flag
                continue
            end

            if path_length > 3
                (path_probs[gs_idxes[gs], gs_idxes[gs_end]]
                 = PathProbs(path_length, dists[gs_end]))
            end
        end
    end
end

function aggregate_paths_probs(path_probs::Matrix{PathProbs}, num_gses:: Int64)
    num_gs_pairs = num_gses*(num_gses-1)/2
    return sum(path.transmissivity for path ∈ path_probs)/(2 * num_gs_pairs)
end

function aggregate_hop_probs(path_probs::Matrix{PathProbs})
    hop_probs = Dict{Int64, HopProb}()
    for path ∈ path_probs
        hop_prob = get!(hop_probs, path.num_nodes, HopProb(0, 0))
        hop_prob.count += 1
        hop_prob.sum_transmissivity += path.transmissivity
    end
    delete!(hop_probs, 0)
    return Dict{Int64, Float64}(
        num => hop_prob.sum_transmissivity / hop_prob.count
        for (num, hop_prob) in hop_probs)
end