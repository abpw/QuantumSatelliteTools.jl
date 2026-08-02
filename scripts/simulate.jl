using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))


using Statistics: mean


using SimpleWeightedGraphs: SimpleWeightedGraph, get_weight
using Graphs: neighbors, vertices
using DataStructures: PriorityQueue

using QuantumSatelliteTools.LossCalculation: total_loss_dB, swapping_loss_dB, dB_to_prob
using QuantumSatelliteTools.Simulator: build_optimized_constellation, Node, Link, simulate,
    NodeLabel, satellite, ground_station
using QuantumSatelliteTools.GenerateGroundStations: generate_grid_gses,
    _generate_population_center_gses

"""
A pair of nodes, for indexing into a path dictionary.
"""
struct NodePair
    node1::Node
    node2::Node

    NodePair(node1, node2) = (
        node1.id < node2.id
        ? new(node1, node2)
        : new(node2, node1)
    )
end

"""
A path between two nodes, consisting of a sequence of links.

# Fields
- `length::Int`: The number of nodes in the path.
- `end_node1::Node`: One endpoint of the path.
- `end_node2::Node`: The other endpoint of the path.
- `inter_nodes::Vector{Node}`: The intermediate nodes in the path, excluding the endpoints.
- `loss_dB::Float64`: The total loss of the path in decibels.
"""
struct Path
    length::Int
    end_node1::Node
    end_node2::Node
    inter_nodes::Vector{Node}
    loss_dB::Float64
end

"""
Generate population ground station nodes with weights.

# Arguments
- `n::Int`: the number of population ground stations to generate.

# Returns
- `Vector{Tuple{GroundStation, Float64}}`: A vector of tuples, each containing a ground
    station (latitude, longitude) and its weight.
"""
function generate_population_center_gs_weights(n::Int)
    rows = _generate_population_center_gses(n)
    gs_pops = [
        (
            (π/180*parse(Float64, row.lat), π/180*parse(Float64, row.lng)),
            parse(Int64, row.population)
        )
        for row ∈ rows
    ]
    max_pop = maximum(pop for (_, pop) ∈ gs_pops)

    return [(gs, pop / max_pop) for (gs, pop) ∈ gs_pops]
end

"""
Custom evaluation function for the grid ground station selection experiment.
Calculates the "importance" of grid ground stations based on their "helpfulness" in
connecting population centers.

# Arguments
- `links::Vector{Link}`: A vector of links to evaluate.

# Returns
- `Dict`: A dictionary mapping each grid ground station's id to its importance score.
"""
function grid_gs_importance_evaluation_fn(links::Vector{Link})
    # Get all unique nodes from the links and index them by their ids
    all_nodes = unique(vcat([link.node1 for link ∈ links], [link.node2 for link ∈ links]))
    node_to_idx = Dict(node.id => i for (i, node) ∈ enumerate(all_nodes))
    idx_to_node = all_nodes

    println("sample link loss: $(total_loss_dB(links[1].channel))")
    return

    # Construct simple weighted graph representation of the network
    sources::Vector{Int} = [node_to_idx[link.node1.id] for link ∈ links]
    destinations::Vector{Int} = [node_to_idx[link.node2.id] for link ∈ links]
    edge_weights::Vector{Float64} = [total_loss_dB(link.channel) for link ∈ links]

    graph::SimpleWeightedGraph{Int,Float64} =
        SimpleWeightedGraph(sources, destinations, edge_weights)

    # Identify population center nodes
    pop_ctr_gses::Vector{Node} = [
        node
        for node ∈ all_nodes
        if node.label.type == ground_station && "aux" ∉ node.label.tags
    ]

    # Construct path matrix to store the paths between population centers
    paths::Dict{NodePair,Path} = Dict{NodePair,Path}()

    # For each population center, perform a modified Dijkstra's to find the
    # paths to other population centers with the highest transmissivity
    for (i, pop_ctr_gs) ∈ enumerate(pop_ctr_gses)
        # Initialize Dijkstra's algorithm data structures
        dists::Vector{Float64} =
            [
                idx_to_node[i] == pop_ctr_gs ? 0 : Inf
                for i ∈ vertices(graph)
            ]
        prevs::Vector{Union{Nothing,Node}} = [nothing for i ∈ vertices(graph)]
        priority_queue::PriorityQueue{Int,Float64} = PriorityQueue{Int,Float64}()

        # Initialize priority queue
        for node_idx ∈ vertices(graph)
            push!(priority_queue, node_idx => dists[node_idx])
        end

        # Iteratively visit the node with the lowest loss
        while !isempty(priority_queue)
            curr_node_idx, curr_dist = popfirst!(priority_queue)
            curr_node = idx_to_node[curr_node_idx]
            node_cost = (
                curr_node ≠ pop_ctr_gs && curr_node.label.type == ground_station
                ? swapping_loss_dB()
                : 0.0
            )

            for neighbor_idx ∈ neighbors(graph, curr_node_idx)
                # println("curr_dist: $(dB_to_prob(curr_dist)), node_cost: $(dB_to_prob(node_cost)), edge_weight: $(dB_to_prob(get_weight(graph, curr_node_idx, neighbor_idx)))")
                new_dist = curr_dist + node_cost + get_weight(graph, curr_node_idx, neighbor_idx)
                # println("New dist: $(dB_to_prob(new_dist))")
                if new_dist < dists[neighbor_idx]
                    dists[neighbor_idx] = new_dist
                    prevs[neighbor_idx] = curr_node
                    priority_queue[neighbor_idx] = new_dist
                end
            end
        end

        # After Dijkstra's, reconstruct the paths to other population centers
        for other_pop_ctr_gs ∈ pop_ctr_gses[(i+1):end]
            # Reconstruct path from pop_ctr_gs to other_pop_ctr_gs
            path_nodes::Vector{Node} = []
            curr_node::Union{Nothing,Node} = other_pop_ctr_gs
            while curr_node ≠ pop_ctr_gs && !isnothing(curr_node)
                push!(path_nodes, curr_node)
                curr_node = prevs[node_to_idx[curr_node.id]]
            end

            # If we reached the source and there are at least 3 nodes
            # (excluding endpoint), we have a valid path
            if curr_node == pop_ctr_gs && length(path_nodes) > 2
                push!(path_nodes, pop_ctr_gs)
                path_loss_dB = dists[node_to_idx[other_pop_ctr_gs.id]]

                paths[NodePair(pop_ctr_gs, other_pop_ctr_gs)] = Path(
                    length(path_nodes),
                    pop_ctr_gs,
                    other_pop_ctr_gs,
                    path_nodes[2:(end-1)],
                    path_loss_dB
                )
            end
        end
    end

    # Use path matrix to calculate importance of each grid ground station
    grid_gs_importance::Dict{Int,Float64} = Dict{Int,Float64}()
    for (pair::NodePair, path::Path) ∈ paths
        # Calculate importance of path
        path_weight = (pair.node1.weight + pair.node2.weight) / 2.0
        path_transmissivity = 1 - dB_to_prob(path.loss_dB)
        println("Path transmissivity: $path_transmissivity")
        path_importance = path_weight * path_transmissivity

        # Assign importance to each grid ground station in the path
        for inter_node::Node ∈ path.inter_nodes
            if (inter_node.label.type == ground_station && "aux" ∈ inter_node.label.tags)
                grid_gs_importance[inter_node.id] = get(grid_gs_importance, inter_node.id, 0.0) + path_importance
            end
        end
    end

    # Return the importance scores for each grid ground station
    return grid_gs_importance

    # Do modified Dijkstra's to find paths between population centers
    # Calculate importance of each path by averaging the weights of the population centers it connects
    # Assign importance to each grid ground station based on the paths it is part of
    # For each aux ground station in a given path, add the path's importance to the station's total importance score
    # Additional considerations:
    # Consider weighting paths by their transmissivity to reflect the quality of the connection they provide
end

"""
Perform the grid ground station selection experiment.

# Arguments
- `duration_s`: total simulation duration in seconds.
- `step_s`: time step in seconds.

# Returns
- Nothing
"""
function main(;
    duration_s=30,
    step_s=30,
)
    # Build a constellation of 1000 satellites based on optimized parameters
    println("Building optimized constellation of 1000 satellites...")
    satellites = build_optimized_constellation()

    # Generate 100 population center ground stations.
    println("Generating 100 population center ground stations...")
    population_center_gses = generate_population_center_gs_weights(100)

    # Generate ground station grids with 200, 400, and 800 km spacing.
    println("Generating ground station grids with 200, 400, and 800 km spacing...")
    grid_gses = Dict(
        200 => generate_grid_gses(equatorial_distance_km=200),
        400 => generate_grid_gses(equatorial_distance_km=400),
        800 => generate_grid_gses(equatorial_distance_km=800),
    )

    # Perform experiments for each grid scenario
    println("Performing experiments for each grid scenario...")
    for (spacing, aux_gses) ∈ grid_gses
        println("Simulating scenario with grid spacing of $spacing km...")

        # Construct the nodes for this scenario
        println("Constructing satellite nodes...")
        sat_nodes = [
            Node(sat, 0.0, NodeLabel(satellite, []))
            for sat ∈ satellites
        ]
        println("Constructing population center ground station nodes...")
        pop_gs_nodes = [
            Node(gs, weight, NodeLabel(ground_station, []))
            for (gs, weight) ∈ population_center_gses
        ]
        println("Constructing grid ground station nodes for this scenario...")
        aux_gs_nodes = [
            Node(gs, 0.0, NodeLabel(ground_station, ["aux"]))
            for gs ∈ aux_gses
        ]
        nodes = vcat(sat_nodes, pop_gs_nodes, aux_gs_nodes)
        println("Total nodes for this scenario: $(length(nodes))")

        # Simulate the scenario and collect metrics
        println("Beginning simulation run for grid spacing $spacing km...")
        metrics = simulate(
            nodes,
            duration_s,
            step_s,
            evaluation_fn=grid_gs_importance_evaluation_fn,
        )

        println("Sample metric: $(metrics[end])")
    end

end