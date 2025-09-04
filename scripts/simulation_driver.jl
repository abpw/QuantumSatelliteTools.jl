using SatelliteToolboxTle
using QuantumSatelliteTools.GenerateGroundStations
using QuantumSatelliteTools.GenerateTLEs
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using QuantumSatelliteTools.AstronomyGeometry: GS
using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight
using Graphs: neighbors, vertices, nv, edges, src, dst
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel
using QuantumSatelliteTools.LossCalculation: reflector_loss, swapping_loss, total_loss, decibels_to_probability
using LinearAlgebra: I
using QuantumSatelliteTools.GenerateGroundStations: generate_population_center_gses
using QuantumSatelliteTools.GenerateTLEs: generate_regular_array_TLEs
using ProgressBars: ProgressBar
using DataStructures: PriorityQueue, enqueue!, dequeue!
using Base: summarysize

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

"""
Helper to add edge data to three edge data arrays.

# Arguments
- `source::Union{OrbitPropagatorSgp4{Float64, Float64},GS}`: A source entity.
- `destination::Union{OrbitPropagatorSgp4{Float64, Float64},GS}`: A destination entity.
- `node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64}`: A node map of entities to integer values.
- `time::Float64`: The time of the transmission.
- `sources::Vector{Int64}`: An array of sources for a graph.
- `destinations::Vector{Int64}`: An array of destinations for a graph.
- `weights::Vector{Float64}`: An array of weights for a graph.

# Returns
- nothing
"""
function edge_data!(source::Union{OrbitPropagatorSgp4{Float64, Float64},GS},
        destination::Union{OrbitPropagatorSgp4{Float64, Float64},GS},
        node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64},
        time::Float64,
        sources::Vector{Int64},
        destinations::Vector{Int64},
        weights::Vector{Float64})
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
                                              representing the satellites in
                                              the network.
- `gses::Vector{GS}`: An array of tuples holding the coordinates of
                                  the ground stations in the network.
- `node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64}`:
                            A map of entities to integer values.
- `time::Float64`: The time that the graph takes place.
- `experiment::Experiment`: The type of experiment.

# Returns
- SimpleWeightedGraph: A graph representing a quantum satellite network with
                       nodes of type Union{OrbitPropagatorSgp4{Float64, Float64},GS}
                       representing satellites and ground stations and edges
                       with weights representing channel transmissivity.
"""
function make_graph(propagators::Vector{OrbitPropagatorSgp4{Float64, Float64}},
        gses::Vector{Tuple{T, T}},
        node_map::Dict{Union{OrbitPropagatorSgp4{Float64, Float64}, GS},Int64},
        time::Float64,
        experiment::Experiment) where T <: Number
    sources, destinations, weights =
        Vector{Int64}(), Vector{Int64}(), Vector{Float64}()
    for propagator ∈ propagators
        for gs ∈ gses
            edge_data!(propagator, gs, node_map, time, sources, destinations, weights)
        end
    end
    if experiment == reflector
        for i ∈ eachindex(propagators)
            for j ∈ i+1:length(propagators)
                edge_data!(propagators[i], propagators[j], node_map, time, sources, destinations, weights)
            end
        end
    end
    return SimpleWeightedGraph(sources, destinations, weights)
end


"""
Conduct a simulation of a quantum satellite network.

# Arguments
- `duration_s::Int64`: The total duration of the simulation in seconds.
- `interval_s::Int64`: The interval at which the system is sampled in seconds.
- `epoch::Float64`: The time at which the simulation begins.
- `gses::Vector{GS}`: An array of ground stations (default GSES).
- `aux_gses::Vector{GS}`: An array of auxiliary groundstations as
                                      produced by the DD experiment
                                      (default missing).
- `experiment::Experiment`: Which experiment to run.
- `num_sats::Int64`: The number of satellites per orbit.

# Returns
- An aggregate representation of the performance of the system under the
  architecture dictated by the experiment.
"""
function simulate(duration_s::Int64, interval_s::Int64, epoch::Float64,
        experiment::Experiment, num_sats::Int64;
        gses::Union{Vector{Tuple{T, T}}, Missing}=missing,
        aux_gses::Union{Vector{Tuple{T, T}}, Missing}=missing,
        time_file::Union{String, Missing}=missing,
        hop_file::Union{String, Missing}=missing,
        agg_file::Union{String, Missing}=missing,
        ) where T <: Number
    # Downloaded Starlink TLEs 8.21.2025
    println("Generating TLEs")
    tles::Vector{TLE} = [tle for θ ∈ 0.0:18.0:162.0 for tle ∈ generate_regular_array_TLEs(orbits=10, sats_per_orbit=num_sats, inclination_rad=θ)]
    propagators = [Propagators.init(Val(:SGP4), tle) for tle ∈ tles]
    gses::Vector{Tuple{Number, Number}} = ismissing(gses) ? generate_population_center_gses(100) : gses
    aux_gses::Vector{Tuple{Number, Number}} = ismissing(aux_gses) ? [] : aux_gses

    # Map entity to an integer index
    println("Mapping nodes to integers")
    node_map = Dict{Union{OrbitPropagatorSgp4{Float64, Float64},GS},Int64}(
        node => index for (index, node) ∈ enumerate(vcat(propagators, gses, aux_gses)))
    # Map integer index of entity to entity type
    println("Mapping node integer values to types")
    types = Dict{Int64, Entities}(
        node_map[propagator] => satellite for propagator ∈ propagators)
    merge!(types, Dict{Int64, Entities}(
        node_map[gs] => ground_station for gs ∈ gses))
    merge!(types, Dict{Int64, Entities}(
        node_map[aux_gs] => aux_ground_station for aux_gs ∈ aux_gses))
    # Perform simulation
    total_prob = 0
    time_transmissivities::Vector{Tuple{Float64, Float64}} = []
    hop_probs = Dict{Int64, HopProb}()
    num_gses = length(gses)
    num_gs_pairs = num_gses*(num_gses-1)/2
    daddy_matrix::Matrix{PathProbs} = [PathProbs(0, 0) for _ in 1:length(gses), _ in 1:length(gses)]
    println("Beginning simulation")
    for time_elapsed ∈ ProgressBar(0:interval_s:duration_s-1)
    # for time_elapsed ∈ 0:interval_s:duration_s-1
        graph = make_graph(propagators, vcat(gses, aux_gses), node_map, epoch+time_elapsed, experiment)
        # println(summarysize(graph))
        all_pairs_path_probs!(
            graph,
            types,
            experiment,
            daddy_matrix
        )
        # println("time elapsed: $(time_elapsed) / 86400")
        # println("hop_probs size: $(summarysize(hop_probs)) bytes")
        # println("graph size: $(summarysize(graph)) bytes")
        # println("path_probs size: $(summarysize(path_probs)) bytes")
        # println("daddy_matrix size: $(summarysize(daddy_matrix)) bytes")
        # curr_hop_probs = aggregate_hop_probs(daddy_matrix)
        for prob ∈ daddy_matrix
            hop_prob = get!(hop_probs, prob.num_nodes, HopProb(0, 0))
            hop_prob.count += 1
            hop_prob.sum_transmissivity += prob.transmissivity
        end
        curr_prob = aggregate_paths_probs(daddy_matrix) / num_gs_pairs
        if !ismissing(time_file)
            push!(time_transmissivities, (time_elapsed, curr_prob))
        end
        total_prob += curr_prob
        if time_elapsed % 3600 == 0 && !ismissing(time_file)
            open(pwd()*"/data/"*time_file, "w") do file
                for transmissivity in time_transmissivities
                    write(file, "$(transmissivity[1]),$(transmissivity[2])\n")
                end
            end
            time_transmissivities = []
            GC.gc()
        end
            for i in 1:length(gses), j in 1:length(gses)
                daddy_matrix[i, j] = PathProbs(0, 0)
            end
    end

    experiment_prob = total_prob / (duration_s ÷ interval_s)
    hop_probs = aggregate_hop_probs(daddy_matrix)
    if !ismissing(hop_file)
        open(pwd()*"/data/"*hop_file, "w") do file
            for hop_prob in hop_probs
                write(file, "$(hop_prob.num_nodes),$(hop_prob.transmissivity)\n")
            end
        end
    end

    if !ismissing(agg_file)
        open(pwd()*"/data/"*agg_file, "w") do file
            write(file, "$(experiment_prob)")
        end
    end
    return experiment_prob
end

function simulation_driver()
    # Perform experiments
    total_probs = Dict{}
    for experiment ∈ [Val(:REF), Val(:DD)]
        simulate
    end

end

function all_pairs_path_probs!(g::SimpleWeightedGraph, types::Dict{Int64, Entities}, experiment::Experiment, path_probs::Matrix{PathProbs})
    gses = [i for i ∈ vertices(g) if types[i] == ground_station]
    gs_idxes = Dict{Int64, Int64}(int => idx for (idx, int) in enumerate(gses))
    node_cost_sat = experiment == reflector ? 1 - reflector_loss() : 1
    node_cost_gs = experiment == reflector ? 1 : 1 - swapping_loss()
    for gs ∈ gses
        priority_queue = PriorityQueue()
        seen = Set{Int64}()
        push!(seen, gs)
        for neighbor ∈ neighbors(g, gs)
            weight = get_weight(g, gs, neighbor)
            enqueue!(priority_queue, Path([gs, neighbor], weight) => -weight)
        end

        while length(priority_queue) > 0
            path = dequeue!(priority_queue)
            curr = path.nodes[end]
            push!(seen, curr)
            if types[curr] == ground_station
                path_probs[gs_idxes[gs], gs_idxes[curr]] = PathProbs(length(path.nodes), path.transmissivity)
            end
            node_cost = types[curr] == satellite ? node_cost_sat : node_cost_gs
            for neighbor in neighbors(g, curr)
                if !(neighbor ∈ seen) && path.transmissivity > 0
                    path_cost = path.transmissivity * node_cost * g.weights[curr, neighbor]
                    enqueue!(priority_queue, Path(vcat(path.nodes, [neighbor]), path_cost) => -path_cost)
                end
            end
        end
    end
    # println(summarysize(path_probs))
    # return aggregate_hop_probs(path_probs)
end

function aggregate_paths_probs(path_probs::Matrix{PathProbs})
    return sum(path.transmissivity for path ∈ path_probs)/2
end

function aggregate_hop_probs(path_probs::Matrix{PathProbs})
    hop_probs = Dict{Int64, Vector{Float64}}()
    for path ∈ path_probs
        hop_prob = get!(hop_probs, path.num_nodes, [0, 0])
        hop_prob[1] += 1
        hop_prob[2] += path.transmissivity
    end
    return [PathProbs(k, s / n) for (k, (s, n)) ∈ hop_probs]
end