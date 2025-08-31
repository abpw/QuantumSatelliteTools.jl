using SatelliteToolboxTle
using QuantumSatelliteTools.GenerateGroundStations
using QuantumSatelliteTools.GenerateTLEs
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using QuantumSatelliteTools.AstronomyGeometry: GS
using SimpleWeightedGraphs, Graphs
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel
using QuantumSatelliteTools.LossCalculation: reflector_loss, swapping_loss, total_loss, decibels_to_probability
using LinearAlgebra: I
using QuantumSatelliteTools.GenerateGroundStations: generate_population_center_gses
using QuantumSatelliteTools.GenerateTLEs: generate_regular_array_TLEs
using ProgressBars: ProgressBar

@enum Entities begin
    satellite
    ground_station
    aux_ground_station
end

@enum Experiment begin
    reflector
    dual_downlink
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
        out_file::Union{String, Missing}=missing) where T <: Number
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
    # println(length(Dict(key => val for (key, val) in types if val == satellite)))
    # println(length(Dict(key => val for (key, val) in types if val == ground_station)))
    # println(length(Dict(key => val for (key, val) in types if val == aux_ground_station)))
    # Perform simulation
    total_prob = 0
    transmissivities::Vector{Tuple{Float64, Float64}} = []
    num_gses = length(gses)
    num_gs_pairs = num_gses*(num_gses-1)/2
    println("Beginning simulation")
    for time_elapsed ∈ ProgressBar(0:interval_s:duration_s-1)
        println("    Time interval: $(time_elapsed)")
        graph = make_graph(propagators, vcat(gses, aux_gses), node_map, epoch+time_elapsed, experiment)
        curr_prob = all_pairs_path_probs(
            graph,
            types,
            experiment) / num_gs_pairs
        println("    Total path probability at interval: $(curr_prob)")
        if !ismissing(out_file)
            push!(transmissivities, (time_elapsed, curr_prob))
        end
        total_prob += curr_prob
    end
    experiment_prob = total_prob / (duration_s ÷ interval_s)
    println("Total probability over experiment: $(experiment_prob)")
    if !ismissing(out_file)
        open(pwd()*"/data/"*out_file, "w") do file
            for transmissivity in transmissivities
                write(file, "$(transmissivity[1]),$(transmissivity[2])\n")
            end
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

function all_pairs_path_probs(g::SimpleWeightedGraph, types::Dict{Int64, Entities}, experiment::Experiment)
    path_probs = zeros(nv(g), nv(g)) + I
    for edge ∈ edges(g)
        source = min(src(edge), dst(edge))
        destination = max(src(edge), dst(edge))
        path_probs[source, destination] = g.weights[source, destination]
    end
    for k ∈ ProgressBar(1:nv(g))
        if experiment == reflector
            if types[k] != satellite
                continue
            end
            node_cost = 1 - reflector_loss()
        else
            node_cost = types[k] == satellite ? 1 : 1 - swapping_loss()
        end
        for i ∈ 1:nv(g)
            if i == k
                continue
            end
            if experiment == dual_downlink && !xor(types[k] == satellite, types[i] == satellite)
                continue
            end
            for j ∈ i+1:nv(g)
                if j == k
                    continue
                end
                if experiment == dual_downlink && !xor(types[k] == satellite, types[j] == satellite)
                    continue
                end
                path_probs[i, j] = max(path_probs[i, j], path_probs[min(i, k), max(i, k)] * path_probs[min(k, j), max(k, j)] * node_cost)
            end
        end
    end
    arr = [path_probs[i, j] for i ∈ vertices(g) if types[i] == ground_station for j ∈ i+1:nv(g) if types[j] == ground_station]
    return sum(arr)
end