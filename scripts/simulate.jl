using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataStructures: PriorityQueue
using Graphs: neighbors, vertices
using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight

using QuantumSatelliteTools.types: Node, Link, satellite, ground_station, NodeLabel, GS
using QuantumSatelliteTools.helpers:
    dB_to_prob, prob_to_dB, _parse_lat_lon, is_ground_station
using QuantumSatelliteTools.GenerateGroundStations:
    generate_population_center_gses, generate_grid_gses
using QuantumSatelliteTools.LossCalculation: swapping_loss_dB, total_loss_dB
using QuantumSatelliteTools.Simulator: build_optimized_constellation, simulate
using QuantumSatelliteTools.VisualizeMap: plot_gs_selections

"""
    NodePair(node1::Node, node2::Node)

A canonically ordered pair of nodes for path indexing.

# Fields

  - `node1::Node`: First node (lower ID).
  - `node2::Node`: Second node (higher ID).
"""
struct NodePair
    node1::Node
    node2::Node

    function NodePair(node1::Node, node2::Node)
        return node1.id < node2.id ? new(node1, node2) : new(node2, node1)
    end
end

"""
    Path

A path between two nodes via a sequence of links.

# Fields

  - `length::Int`: Number of nodes in the path.
  - `end_node1::Node`: First endpoint.
  - `end_node2::Node`: Second endpoint.
  - `inter_nodes::Vector{Node}`: Intermediate nodes (excluding endpoints).
  - `loss_dB::Float64`: Total path loss [dB].
"""
struct Path
    length::Int
    end_node1::Node
    end_node2::Node
    inter_nodes::Vector{Node}
    loss_dB::Float64
end

"""
    is_aux_gs(node::Node) -> Bool

Check if a `node` is an auxiliary ground station.
"""
function is_aux_gs(node::Node)
    return is_ground_station(node) && :aux ∈ node.label.tags
end

"""
    is_non_aux_gs(node::Node) -> Bool

Check if a `node` is a non-auxiliary ground station.
"""
function is_non_aux_gs(node::Node)
    return is_ground_station(node) && :aux ∉ node.label.tags
end

"""
    generate_population_center_gs_weights(n::Integer) -> Vector{Tuple{GS, Float64}}

Generate `n` population center ground stations with normalized weights.

# Returns

  - `Vector{Tuple{GS, Float64}}`: Ground stations and their normalized population weights.
"""
function generate_population_center_gs_weights(n::Integer)
    rows = generate_population_center_gses(Val(:rows), n)
    gs_populations = [
        (_parse_lat_lon(row.lat, row.lng), parse(Int64, row.population)) for row in rows
    ]
    max_population = maximum(population for (_, population) in gs_populations)

    return [(gs, population / max_population) for (gs, population) in gs_populations]
end

"""
    grid_gs_importance_evaluation_fn(links::Vector{Link}) -> Dict{Node, Real}

Evaluate grid ground station importance in connecting population centers for the network
state given by `links`.

Computes importance scores based on path quality and frequency of usage as intermediate
nodes.

# Returns

  - `Dict{Node, Real}`: Importance scores for grid ground stations.
"""
function grid_gs_importance_evaluation_fn(links::Vector{Link})
    all_nodes = unique(vcat([link.node1 for link in links], [link.node2 for link in links]))
    node_to_idx::Dict{Int,Int} = Dict(node.id => i for (i, node) in enumerate(all_nodes))
    idx_to_node::Vector{Node} = all_nodes

    sources::Vector{Int} = [node_to_idx[link.node1.id] for link in links]
    destinations::Vector{Int} = [node_to_idx[link.node2.id] for link in links]
    weights::Vector{Float64} = [total_loss_dB(link.channel) for link in links]

    graph = SimpleWeightedGraph(sources, destinations, weights)

    pop_center_gses::Vector{Node} = [node for node in all_nodes if is_non_aux_gs(node)]

    paths = Dict{NodePair,Path}()

    for (i, pop_center_gs) in enumerate(pop_center_gses)
        dists::Vector{Float64} = [
            idx_to_node[i] == pop_center_gs ? 0.0 : Inf for i in vertices(graph)
        ]
        prevs::Vector{Union{Nothing,Int}} = [nothing for i in vertices(graph)]
        pq = PriorityQueue{Int,Float64}()

        for node_idx in vertices(graph)
            push!(pq, node_idx => dists[node_idx])
        end

        while !isempty(pq)
            curr_node_idx, curr_dist = popfirst!(pq)
            curr_node = idx_to_node[curr_node_idx]

            is_swap_gs = curr_node ≠ pop_center_gs && is_ground_station(curr_node)
            node_cost = is_swap_gs ? swapping_loss_dB() : 0.0

            for neighbor_idx in neighbors(graph, curr_node_idx)
                new_dist =
                    curr_dist + node_cost + get_weight(graph, curr_node_idx, neighbor_idx)

                if new_dist < dists[neighbor_idx]
                    dists[neighbor_idx] = new_dist
                    prevs[neighbor_idx] = curr_node_idx
                    pq[neighbor_idx] = new_dist
                end
            end
        end

        for other_pop_center_gs in pop_center_gses[(i+1):end]
            path_nodes::Vector{Node} = []
            curr_idx = node_to_idx[other_pop_center_gs.id]
            source_idx = node_to_idx[pop_center_gs.id]
            while curr_idx ≠ source_idx && !isnothing(prevs[curr_idx])
                push!(path_nodes, idx_to_node[curr_idx])
                curr_idx = prevs[curr_idx]
            end

            if curr_idx == source_idx && length(path_nodes) > 2
                push!(path_nodes, pop_center_gs)
                path_loss_dB = dists[node_to_idx[other_pop_center_gs.id]]

                paths[NodePair(pop_center_gs, other_pop_center_gs)] = Path(
                    length(path_nodes),
                    pop_center_gs,
                    other_pop_center_gs,
                    path_nodes[2:(end-1)],
                    path_loss_dB,
                )
            end
        end
    end

    grid_gs_importance = Dict{Node,Float64}()
    for (pair::NodePair, path::Path) in paths
        path_weight = (pair.node1.weight + pair.node2.weight) / 2.0
        path_transmissivity = 1 - dB_to_prob(path.loss_dB)
        path_importance = path_weight * path_transmissivity

        for inter_node in path.inter_nodes
            if is_aux_gs(inter_node)
                grid_gs_importance[inter_node] =
                    get(grid_gs_importance, inter_node, 0.0) + path_importance
            end
        end
    end

    return grid_gs_importance
end

"""
    visualize_gs_importances(
        spacing::Integer, pop_center_gses::Vector{Tuple{GS,Real}},
        metrics::Vector{Dict{Node,Real}}
    ) -> Nothing

Visualize grid ground station importance on a world map.

# Arguments

  - `spacing::Integer`: Grid spacing [km].
  - `pop_center_gses::Vector{Tuple{GS,Float64}}`: Population centers with weights.
  - `metrics::Vector{Dict{Node,Real}}`: Importance metrics by time step.
"""
function visualize_gs_importances(
    spacing::Integer,
    pop_center_gses::Vector{Tuple{GS,Float64}},
    metrics::Vector{Dict{Node,T}},
) where {T<:Real}
    importances = Dict{GS,Float64}()
    for time_step in metrics
        for (node, importance) in time_step
            importances[node.obj] = get(importances, node.obj, 0.0) + importance
        end
    end

    filtered_importances = Dict(
        gs => importance for (gs, importance) in importances if importance > 0.0
    )

    plot_gs_selections(spacing, pop_center_gses, filtered_importances)
end

"""
    main(; duration_s::Real=86400, step_s::Real=300) -> Nothing

Execute grid ground station selection simulation experiment.

Build a 1000-satellite constellation, evaluate auxiliary ground station importance
across three grid spacings (200/400/800 km), and visualize results.

Arguments

  - `duration_s::Real=86400`: Simulation duration [s].
  - `step_s::Real=300`: Time step [s].
"""
function main(; duration_s::Real=86400, step_s::Real=300)
    println("Building optimized constellation of 1000 satellites...")
    satellites = build_optimized_constellation()

    println("Generating 100 population center ground stations...")
    pop_center_gses::Vector{Tuple{GS,Float64}} = generate_population_center_gs_weights(100)

    println("Generating ground station grids with 200, 400, and 800 km spacing...")
    grid_gses = Dict(
        200 => generate_grid_gses(; equatorial_distance_km=200, force_land=true),
        400 => generate_grid_gses(; equatorial_distance_km=400, force_land=true),
        800 => generate_grid_gses(; equatorial_distance_km=800, force_land=true),
    )

    println("Performing experiments for each grid scenario...")
    for (spacing::Int, aux_gses::Vector{GS}) in grid_gses
        println("Simulating scenario with grid spacing of $spacing km...")

        println("Constructing satellite nodes...")
        sat_nodes = [Node(sat, 0, NodeLabel(satellite)) for sat in satellites]

        println("Constructing population center ground station nodes...")
        pop_gs_nodes = [
            Node(gs, weight, NodeLabel(ground_station)) for (gs, weight) in pop_center_gses
        ]

        println("Constructing grid ground station nodes for this scenario...")
        aux_gs_nodes = [
            Node(aux_gs, 0, NodeLabel(ground_station, Set([:aux]))) for aux_gs in aux_gses
        ]

        nodes = vcat(sat_nodes, pop_gs_nodes, aux_gs_nodes)
        println("Total nodes for this scenario: $(length(nodes))")

        println("Beginning simulation run for grid spacing $spacing km...")
        metrics::Vector{Dict{Node,Float64}} = simulate(
            nodes, duration_s, step_s; evaluation_fn=grid_gs_importance_evaluation_fn
        )

        println("Visualizing results for grid spacing $spacing km...")
        visualize_gs_importances(spacing, pop_center_gses, metrics)
    end
end