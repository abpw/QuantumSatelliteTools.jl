using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using ..GenerateSatellites: generate_regular_array_TLEs

@enum NodeLabels begin
    sat
    gs
    aux_gs
end

struct Node
    id::Int
    obj::Union{Propagator, GroundStation}
    weight::Float64
    labels::Vector{NodeLabels}
end

# TODO: write Channel to abstract Freespace/Fiber stuff

struct Link
    id::Int
    node1::Node
    node2::Node
    freespace_channel::FreespaceChannel
end

# TODO: write pairwise function for downlink

"""
Build a synthetic constellation and initialize SGP4 propagators based on the
optimized inclination angles and satellite allocations discovered in
https://arxiv.org/abs/2603.02480. The five base orbits have inclination angles 
of 153.21, 24.72, 143.04, 69.94, and 149.05 degrees with 45, 37, 15, 1, and 2
satellites respectively.

# Keyword Arguments
- `orbits`: Number of orbits per inclination angle (default 10).

# Returns
- A vector of SGP4 propagators initialized with the generated TLEs.
"""
function build_optimized_constellation(;
    orbits::Int=10,
)
    satellite_parameters = zip(
        deg2rad.([153.21, 24.72, 143.04, 69.94, 149.05]),
        [45, 37, 15, 1, 2],
    )

    propagators = Propagator[]
    
    for (inclination_rad, sats_per_orbit) in satellite_parameters
        tles = generate_regular_array_TLEs(
            orbits=orbits,
            sats_per_orbit=sats_per_orbit,
            altitude_km=550,
            inclination_rad=inclination_rad,
        )

        append!(
            propagators,
            [Propagators.init(Val(:SGP4), tle) for tle in tles]
        )  
    end
    return propagators
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
function simulate_step(;
    nodes::Vector{Node},
    PairwiseEvaluationFn::Function,
    evaluation_fn::Function,
    time::Float64,
    tx_aperture_m::Float64=0.6,
    rx_aperture_m::Float64=0.6,
    wavelength_nm::Float64=1550.0,
)
    # get freespace channels for all satellite-ground station pairs
    channels = FreespaceChannel[]
    for node in nodes
        if NodeLabels.sat in node.labels
            for gs_node in nodes
                if NodeLabels.gs in gs_node.labels
                    channel = FreespaceChannel(node.obj, gs_node.obj, transmitter_diameter_m=tx_aperture_m, receiver_diameter_m=rx_aperture_m, wavelength_nm=wavelength_nm, time=time)
                    if !isnothing(channel)
                        push!(channels, channel)
                    end
                end
            end
        end
    end

    # apply the evaluation function to the current state of the channels
    return evaluation_fn(channels)
end

"""
Simulate the constellation performance over a specified duration and time step.
"""
function simulate(
    propagators::Vector{Propagator},
    gses::Vector{GroundStation},
    aux_gses::Vector{GroundStation}=[],
    duration_s::Int,
    step_s::Int;
    tx_aperture_m::Float64=0.6,
    rx_aperture_m::Float64=0.6,
    wavelength_nm::Float64=1550.0,
)
    for t in 0:step_s:duration_s
        # build nodes for this time step
        # nodes = Node[]

        # for (i, sat) in enumerate(propagators)
        #     push!(nodes, Node(i, sat, 0, [NodeLabels.sat]))
        # end

        # for (i, gs) in enumerate(gses)
        #     # TODO: assign weights based on population
        #     push!(nodes, Node(length(propagators) + i, gs, 1.0, [NodeLabels.gs]))
        # end

        # for (i, aux_gs) in enumerate(aux_gses)
        #     push!(nodes, Node(length(propagators) + length(gses) + i, aux_gs, 0, [NodeLabels.aux_gs]))
        # end

        # evaluate the current state of the constellation
        metrics = simulate_step(
            nodes=nodes,
            evaluation_fn=evaluate_constellation,
            time=t,
            tx_aperture_m=tx_aperture_m,
            rx_aperture_m=rx_aperture_m,
            wavelength_nm=wavelength_nm,
        )

    end
end
