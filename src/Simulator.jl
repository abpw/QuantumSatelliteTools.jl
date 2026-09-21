using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using SatelliteToolboxBase: JD_J2000

using ..constants: SECONDS_PER_DAY
using ..types: GS, Propagator, Node, Link
using ..helpers: is_satellite, is_ground_station, _seconds_to_JD
using ..helpers: validate_non_negative_finite
using ..FreespaceChannels: AbstractChannel, FreespaceChannel
using ..GenerateSatellites: generate_regular_array_TLEs
using ..LossCalculation: total_loss

"""
    build_optimized_constellation(; orbits::Integer=10) -> Vector{Propagator}

Build synthetic constellation and initialize SGP4 propagators based on
optimized inclination angles and satellite allocations discovered in
https://arxiv.org/abs/2603.02480.

The five base orbits have inclination angles of 153.21, 24.72, 143.04, 69.94,
and 149.05 degrees with 45, 37, 15, 1, and 2 satellites respectively.

# Keyword Arguments

  - `orbits::Integer=10`: Number of orbits per inclination angle.

# Returns

  - `Vector{Propagator}`: SGP4 propagators initialized with generated TLEs.
"""
function build_optimized_constellation(; orbits::Integer=10)
    satellite_parameters = zip(
        deg2rad.([153.21, 24.72, 143.04, 69.94, 149.05]), [45, 37, 15, 1, 2]
    )

    propagators = Propagator[]

    for (inclination_rad, sats_per_orbit) in satellite_parameters
        tles = generate_regular_array_TLEs(;
            orbital_planes=orbits,
            sats_per_plane=sats_per_orbit,
            altitude_km=1000,
            inclination_rad,
        )

        append!(propagators, [Propagators.init(Val(:SGP4), tle) for tle in tles])
    end
    return propagators
end

"""
    create_link(node1, node2, pairwise_fns; kwargs...) -> Union{Link, Nothing}

Create link between two nodes if pairwise conditions are met.

# Arguments

  - `node1::Node`: First node.
  - `node2::Node`: Second node.
  - `pairwise_fns::Vector{<:Function}`: Functions that take two nodes and return true if
    conditions are met, or false otherwise.

# Keyword Arguments

  - `kwargs...`: Additional keyword arguments to pass to Link.

# Returns

  - `Union{Link, Nothing}`: New link if conditions are met, or `nothing` otherwise.
"""
function create_link(node1::Node, node2::Node, pairwise_fns::Vector{<:Function}; kwargs...)
    for pairwise_fn in pairwise_fns
        if pairwise_fn(node1, node2)
            channel = FreespaceChannel(
                node1.obj, node2.obj; kwargs...
            )
            if (!isnothing(channel))
                return Link(node1, node2, channel)
            end
        end
    end
    return nothing
end

"""
    downlink_pairwise_fn(node1::Node, node2::Node) -> Bool

Check if `node1` and `node2` are a satellite and a ground station.
"""
function downlink_pairwise_fn(node1::Node, node2::Node)
    return (is_satellite(node1) && is_ground_station(node2)) ||
        (is_ground_station(node1) && is_satellite(node2))
end

"""
    intersat_pairwise_fn(node1::Node, node2::Node) -> Bool

Check if `node1` and `node2` are both satellites.
"""
function intersat_pairwise_fn(node1::Node, node2::Node)
    return is_satellite(node1) && is_satellite(node2)
end

"""
    default_evaluation_fn(links::Vector{Link}) -> Float64

Compute total transmissivity across all `links`.
"""
function default_evaluation_fn(links::Vector{Link})
    if isempty(links)
        return 0.0
    end
    return sum(1 - total_loss(link.channel) for link in links)
end

"""
    simulate_step(
        pairwise_fns::Vector{<:Function}, evaluation_fn::Function,
        nodes::Vector{Node}, time::Real;
        transmitter_diameter_m::Real=0.6, receiver_diameter_m::Real=0.6,
        wavelength_nm::Real=1550.0,
    ) -> Any

Perform a single time step of the simulation, applying the evaluation function.

# Arguments

  - `pairwise_fns::Vector{<:Function}`: Functions with signature `(::Node, ::Node) -> Bool`
    to identify valid pairs of nodes for link creation.
  - `evaluation_fn`: Function with signature `(::Vector{Link}) -> Any` to evaluate current
    state of network links.
  - `nodes::Vector{Node}`: Nodes in network (satellites and ground stations).
  - `time::Real`: Current time step [s].

# Keyword Arguments

  - `transmitter_diameter_m::Real=0.6`: Transmitting telescope diameter [m].
  - `receiver_diameter_m::Real=0.6`: Receiving telescope diameter [m].
  - `wavelength_nm::Real=1550.0`: Wavelength [nm].

# Returns

  - `Any`: The result of the evaluation function.
"""
function simulate_step(
    pairwise_fns::Vector{<:Function},
    evaluation_fn::Function,
    nodes::Vector{Node},
    time::Real;
    transmitter_diameter_m::Real=0.6,
    receiver_diameter_m::Real=0.6,
    wavelength_nm::Real=1550.0,
)
    links::Vector{Link} = []
    for i in 1:length(nodes)
        for j in (i+1):length(nodes)
            link = create_link(
                nodes[i],
                nodes[j],
                pairwise_fns;
                time,
                transmitter_diameter_m,
                receiver_diameter_m,
                wavelength_nm,
            )
            if !isnothing(link)
                push!(links, link)
            end
        end
    end

    return evaluation_fn(links)
end

"""
    simulate(nodes::Vector{Node}, duration_s::Integer, step_s::Integer) -> Vector{Any}

Simulate network behavior over a specified duration and time step.

# Arguments

  - `nodes::Vector{Node}`: Nodes in network (satellites and ground stations).
  - `duration_s::Integer`: Total simulation duration [s].
  - `step_s::Integer`: Time step for simulation [s].

# Keyword Arguments

  - `epoch::Real=JD_J2000`: Epoch time for simulation [Julian date].
  - `pairwise_fns::Vector{<:Function}=[downlink_pairwise_fn]`: Functions with signature
    `(::Node, ::Node) -> Bool` to identify valid pairs of nodes for link creation.
  - `evaluation_fn::Function=default_evaluation_fn`: Function with signature
    `(::Vector{Link}) -> Any` to evaluate current state of network links.
  - `transmitter_diameter_m::Real=0.6`: Transmitting telescope diameter [m].
  - `receiver_diameter_m::Real=0.6`: Receiving telescope diameter [m].
  - `wavelength_nm::Real=1550.0`: Wavelength [nm].

# Returns

  - `Vector{Any}`: Results of the evaluation function at each time step.

# Throws

  - `ArgumentError`: If any of the input parameters are invalid.
"""
function simulate(
    nodes::Vector{Node},
    duration_s::Integer,
    step_s::Integer;
    epoch::Real=JD_J2000,
    pairwise_fns::Vector{<:Function}=[downlink_pairwise_fn],
    evaluation_fn::Function=default_evaluation_fn,
    transmitter_diameter_m::Real=0.6,
    receiver_diameter_m::Real=0.6,
    wavelength_nm::Real=1550.0,
)
    validate_non_negative_finite("duration_s", duration_s)
    validate_non_negative_finite("step_s", step_s)
    validate_non_negative_finite("transmitter_diameter_m", transmitter_diameter_m)
    validate_non_negative_finite("receiver_diameter_m", receiver_diameter_m)
    validate_non_negative_finite("wavelength_nm", wavelength_nm)
    if !isfinite(epoch)
        throw(ArgumentError("epoch must be finite, got $epoch"))
    end

    metrics::Vector{Any} = []
    for t in 0.0:step_s:duration_s
        println("Simulating time step $t seconds...")
        push!(
            metrics,
            simulate_step(
                pairwise_fns,
                evaluation_fn,
                nodes,
                epoch + _seconds_to_JD(t);
                transmitter_diameter_m,
                receiver_diameter_m,
                wavelength_nm,
            ),
        )
    end

    return metrics
end
