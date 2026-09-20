using Dates: DateTime

using StaticArrays: SVector
using SatelliteToolboxPropagators: OrbitPropagatorSgp4

############################################################################################
#                                     Coordinate Types                                     #
############################################################################################

"""
    InputGS

Ground station coordinates with any Real precision.
"""

const InputGS = Union{
    Tuple{<:Real,<:Real},
    Tuple{<:Real,<:Real,<:Real},
    SVector{2,<:Real},
    SVector{3,<:Real}
}

"""
    GS

Ground station coordinates with Float64 precision.
"""

const GS = Union{
    NTuple{2,Float64},
    NTuple{3,Float64},
    SVector{2,Float64},
    SVector{3,Float64}
}

"""
    Point3D

A 3D point.
"""
const Point3D = Union{
    Tuple{<:Real,<:Real,<:Real},
    SVector{3,<:Real}
}

############################################################################################
#                                      Network Types                                       #
############################################################################################

"""
    AbstractChannel

Abstract base type for communication channels.
"""
abstract type AbstractChannel end

"""
    Conditions

Atmospheric conditions affecting signal propagation.

# Variants

  - `clear`: Clear sky (value: 0).
  - `fog`: Fog conditions (value: 1).
  - `rain`: Rain conditions (value: 2).
  - `snow`: Snow conditions (value: 3).
"""
@enum Conditions begin
    clear
    fog
    rain
    snow
end

"""
    LightCondition

Lighting condition of a point relative to the Sun and Earth's shadow.

# Variants
  - `umbra`: Fully in Earth's shadow; no direct sunlight (value: 0).
  - `penumbra`: Partially in Earth's shadow (value: 1).
  - `sunlight`: Fully illuminated by the Sun (value: 2).
"""
@enum LightCondition begin
    umbra
    penumbra
    sunlight
end

"""
    Propagator
    
Orbit propagator using SGP4 with Float64 precision.
"""
const Propagator = OrbitPropagatorSgp4{Float64,Float64}

"""
    NodeTypes

Roles or types of a node in the quantum network.

# Variants

  - `satellite`: Satellite node (value: 0).
  - `ground_station`: Ground station node (value: 1).
"""
@enum NodeTypes begin
    satellite
    ground_station
end

"""
    NodeLabel

Label for a node in the network graph.

# Fields

  - `type::NodeTypes`: Type of the node (satellite or ground station).
  - `tags::Set{Symbol}`: Optional tags for the node (e.g., :aux).

# Constructors

    NodeLabel(type::NodeTypes, tags::Set{Symbol}=Set{Symbol}()) -> NodeLabel

Create a node label with the given `type` and optional `tags`.
"""
struct NodeLabel
    type::NodeTypes
    tags::Set{Symbol}

    NodeLabel(type::NodeTypes, tags::Set{Symbol}=Set{Symbol}()) = new(type, tags)
end

"""
    _make_id_counter() -> Function

Create an ID counter closure that generates unique sequential integers starting from 1.
"""
function _make_id_counter()
    count = 0
    return () -> (count += 1)
end

"""
    _next_node_id()

Get the next unique node ID.
"""
_next_node_id = _make_id_counter()

"""
    _next_link_id()

Get the next unique link ID.
"""
_next_link_id = _make_id_counter()

"""
    Node

Node in network graph representing either a propagator or ground station.

# Fields

  - `id::Int`: Unique identifier (auto-generated).
  - `obj::Union{Propagator, GS}`: Underlying satellite or ground station object.
  - `weight::Float64`: Weight associated with the node (used in graph algorithms).
  - `label::NodeLabel`: Label describing the node's properties.

# Constructors

    Node(obj::Union{Propagator, GS}, weight::Real, label::NodeLabel) -> Node

Create a node with the given `obj`, `weight`, and `label`. The node ID is automatically
generated.
"""
struct Node
    id::Int
    obj::Union{Propagator,GS}
    weight::Float64
    label::NodeLabel

    function Node(obj::Union{Propagator,GS}, weight::Real, label::NodeLabel)
        new(_next_node_id(), obj, Float64(weight), label)
    end
end

"""
    Link{C<:AbstractChannel}

Communication link between two nodes over a free-space channel.

# Type Parameters

  - `C<:AbstractChannel`: Type of the channel used for communication.

# Fields

  - `id::Int`: Unique identifier (auto-generated).
  - `node1::Node`: First endpoint of the link.
  - `node2::Node`: Second endpoint of the link.
  - `channel::C`: Model of the link's physical channel.

# Constructors

    Link(node1::Node, node2::Node, channel::C) where {C<:AbstractChannel} -> Link

Create a link between `node1` and `node2` with the specified `channel`. The link ID is
automatically generated.
"""
struct Link{C<:AbstractChannel}
    id::Int
    node1::Node
    node2::Node
    channel::C

    function Link(node1::Node, node2::Node, channel::C) where {C<:AbstractChannel}
        return new{C}(_next_link_id(), node1, node2, channel)
    end
end

############################################################################################
#                                       Other Types                                        #
############################################################################################

"""
    Time

Time representation as either Julian days (Real), DateTime, or Missing.
"""
const Time = Union{<:Real,DateTime,Missing}