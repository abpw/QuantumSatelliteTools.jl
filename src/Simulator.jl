using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using ..GenerateSatellites: generate_regular_array_TLEs
using ..FreespaceChannels: FreespaceChannel, AbstractChannel
using ..LossCalculation: total_loss

# Ground stations are represented as (latitude_rad, longitude_rad) tuples.
const GroundStation = Tuple{Float64,Float64}
const Propagator = OrbitPropagatorSgp4{Float64,Float64}

"""
Enum identifying the role or type of a node in the quantum network.

# Variants
- `sat`: A satellite node.
- `gs`: A primary ground station node.
- `aux_gs`: An auxiliary ground station node.
"""
@enum NodeTypes begin
    satellite
    ground_station
end

"""
A label for a node in the network graph.

# Fields
- `type::NodeTypes`: The type of the node (e.g., satellite, ground station).
- `tags::Vector{String}`: Additional tags for the node
    (e.g., aux for ground stations).
"""
struct NodeLabel
    type::NodeTypes
    tags::Vector{String}
end

"""
Function for making an ID counter closure.
"""
function make_id_counter()
    count = 0
    return () -> (count += 1)
end

next_node_id = make_id_counter()
next_link_id = make_id_counter()

"""
A node in the network graph representing either a propagator or ground station.

# Fields
- `id::Int`: Unique identifier for the node.
- `obj::Union{Propagator, GroundStation}`: The underlying object.
- `weight::Float64`: Weight associated with the node, used in graph algorithms.
- `labels::Vector{NodeLabels}`: Collection of labels describing node properties.
"""
struct Node
    id::Int
    obj::Union{Propagator,GroundStation}
    weight::Float64
    label::NodeLabel

    Node(obj, weight, label) = new(next_node_id(), obj, weight, label)
end

"""
A communication link between two nodes over a freespace channel.

# Fields
- `id::Int`: Unique identifier for the link.
- `node1::Node`: One endpoint of the link.
- `node2::Node`: The other endpoint of the link.
- `channel:: AbstractChannel`: The model of the link's physical channel.
"""
struct Link
    id::Int
    node1::Node
    node2::Node
    channel::AbstractChannel

    Link(node1, node2, channel) = new(next_link_id(), node1, node2, channel)
end

"""
Build a synthetic constellation and initialize SGP4 propagators based on the
optimized inclination angles and satellite allocations discovered in
https://arxiv.org/abs/2603.02480. The five base orbits have inclination angles 
of 153.21, 24.72, 143.04, 69.94, and 149.05 degrees with 45, 37, 15, 1, and 2
satellites respectively.

# Keyword Arguments
- `orbits`: Number of orbits per inclination angle (default 10).

# Returns
- `Vector{Propagator}`: A vector of SGP4 propagators initialized with the
    generated TLEs.
"""
function build_optimized_constellation(;
    orbits::Int=10,
)
    satellite_parameters = zip(
        deg2rad.([153.21, 24.72, 143.04, 69.94, 149.05]),
        [45, 37, 15, 1, 2],
    )

    propagators = Propagator[]

    for (inclination_rad, sats_per_orbit) ∈ satellite_parameters
        tles = generate_regular_array_TLEs(
            orbital_planes=orbits,
            sats_per_plane=sats_per_orbit,
            altitude_km=1000,
            inclination_rad=inclination_rad,
        )

        append!(
            propagators,
            [Propagators.init(Val(:SGP4), tle) for tle ∈ tles]
        )
    end
    return propagators
end

"""
Create a link between two nodes if the pairwise conditions are met.

# Arguments
- `node1::Node`: First node.
- `node2::Node`: Second node.
- `pairwise_fn::Function`: A function that takes two nodes and returns true if
    conditions are met, or false otherwise.

# Keyword Arguments
- `kwargs...`: Additional keyword arguments to pass to Link.

# Returns
- `Union{Link, Nothing}`: The created link if conditions are met,
    or `nothing` otherwise.
"""
function create_link(node1::Node, node2::Node, pairwise_fns::Vector{<:Function}; kwargs...)
    for pairwise_fn ∈ pairwise_fns
        if pairwise_fn(node1, node2)
            channel::Union{FreespaceChannel,Nothing} = FreespaceChannel(node1.obj, node2.obj; kwargs...)
            if (!isnothing(channel))
                return Link(node1, node2, channel)
            end
        end
    end
    return nothing
end

"""
Dual downlink pairwise function that checks if one node is a satellite and the
other is a ground station.

# Arguments
- `node1::Node`: First node.
- `node2::Node`: Second node.

# Returns
- `Bool`: `true` if one node is a satellite and the other is a ground station,
  `false` otherwise.
"""
function downlink_pairwise_fn(node1::Node, node2::Node)
    return ((node1.label.type == satellite && node2.label.type == ground_station) ||
            (node1.label.type == ground_station && node2.label.type == satellite))
end

"""
Inter-satellite pairwise function that checks if both nodes are satellites.

# Arguments
- `node1::Node`: First node.
- `node2::Node`: Second node.

# Returns
- `Bool`: `true` if both nodes are satellites, `false` otherwise.
"""
function intersat_pairwise_fn(node1::Node, node2::Node)
    return node1.label.type == satellite && node2.label.type == satellite
end

"""
Default evaluation function that computes the total transmissivity across all
links.

# Arguments
- `links::Vector{Link}`: A vector of links to evaluate.

# Returns
- `Float64`: The total transmissivity across all links.
"""
function default_evaluation_fn(links::Vector{Link})
    if isempty(links)
        return 0.0
    end
    return sum(1 - total_loss(link.channel) for link ∈ links)
end

"""
Perform a single time step of the simulation, applying the evaluation function.

# Keyword Arguments
- `nodes`: vector of nodes (satellites and ground stations).
- `evaluation_fn`: function that takes the current set of freespace channels
    and their endpoints and returns a metric.
- `time`: current time step in seconds.

# Returns
- The result of the evaluation function.
"""
function simulate_step(
    pairwise_fns::Vector{<:Function},
    evaluation_fn::Function,
    nodes::Vector{Node},
    time::Float64;
    tx_aperture_m::Float64=0.6,
    rx_aperture_m::Float64=0.6,
    wavelength_nm::Float64=1550.0,
)
    # get links for all valid pairs of nodes
    links = Link[]
    for i ∈ 1:length(nodes)
        for j ∈ (i+1):length(nodes)
            link = create_link(
                nodes[i],
                nodes[j],
                pairwise_fns,
                time=time,
                transmitter_diameter_m=tx_aperture_m,
                receiver_diameter_m=rx_aperture_m,
                wavelength_nm=wavelength_nm
            )
            if !isnothing(link)
                push!(links, link)
            end
        end
    end
    # apply the evaluation function to the current state of the channels
    return evaluation_fn(links)
end

"""
Simulate the constellation performance over a specified duration and time step.
"""
function simulate(
    nodes::Vector{Node},
    duration_s::Int,
    step_s::Int;
    pairwise_fns::Vector{<:Function}=[downlink_pairwise_fn],
    evaluation_fn::Function=default_evaluation_fn,
    tx_aperture_m::Float64=0.6,
    rx_aperture_m::Float64=0.6,
    wavelength_nm::Float64=1550.0,
)
    # collect metrics at each time step
    metrics = []

    for t ∈ 0.0:step_s:duration_s
        # evaluate the current state of the network
        println("Simulating time step $t seconds...")
        push!(metrics, simulate_step(
            pairwise_fns,
            evaluation_fn,
            nodes,
            t,
            tx_aperture_m=tx_aperture_m,
            rx_aperture_m=rx_aperture_m,
            wavelength_nm=wavelength_nm,
        ))
    end

    return metrics
end
