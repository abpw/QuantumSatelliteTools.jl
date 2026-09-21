using SatelliteAnalysis: is_ground_facility_visible
using SatelliteToolboxBase: OrbitStateVector
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using SatelliteToolboxTransformations: ecef_to_ned

using ..constants: SEMIMAJOR_RADIUS, SEMIMINOR_RADIUS
using ..helpers: norm, _get_time_JD, _JD_to_seconds
using ..types: AbstractChannel, GS, Time, LightCondition, umbra, Conditions, clear
using ..AstronomyGeometry: light_at_point
using ..AstronomyGeometry: ellipsoid_line_intersection, geodetic_to_ecef, eci_to_ecef

"""
    FreespaceChannel <: AbstractChannel

Free-space optical channel with atmospheric and geometric properties.

# Fields

  - `distance_m::Float64`: Channel distance [m].
  - `elevation_angle_rad::Float64`: Elevation angle [rad].
  - `min_altitude_m::Float64`: Minimum altitude [m].
  - `transmitter_diameter_m::Float64`: Transmitting telescope diameter [m].
  - `receiver_diameter_m::Float64`: Receiving telescope diameter [m].
  - `wavelength_nm::Float64`: Wavelength [nm].
  - `conditions::Conditions`: Atmospheric conditions.
  - `light_condition::LightCondition`: Lighting condition.

# Constructors

    FreespaceChannel(
        distance_m::Real, elevation_angle_rad::Real;
        min_altitude_m::Real=0, transmitter_diameter_m::Real=0.6,
        receiver_diameter_m::Real=0.6, wavelength_nm::Real=1550,
        conditions::Conditions=clear, light_condition::LightCondition=umbra,
    ) -> FreespaceChannel
    FreespaceChannel(
        sat::OrbitPropagatorSgp4, gs::GS;
        transmitter_diameter_m::Real=0.6, receiver_diameter_m::Real=0.6,
        wavelength_nm::Real=1550, time::Time=missing,
        conditions::Conditions=clear, min_angle::Real=deg2rad(20),
    ) -> Union{FreespaceChannel, Nothing}
    FreespaceChannel(
        sat1::OrbitPropagatorSgp4, sat2::OrbitPropagatorSgp4;
        transmitter_diameter_m::Real=0.6, receiver_diameter_m::Real=0.6,
        wavelength_nm::Real=1550, time::Time=missing, conditions::Conditions=clear,
    ) -> Union{FreespaceChannel, Nothing}

Create a `FreespaceChannel` if visibility constraint is met.

The first constructor creates a channel with specified parameters.
The second constructor creates a channel between a satellite and a ground station.
The third constructor creates a channel between two satellites.

# Arguments

  - `sat::OrbitPropagatorSgp4`: Satellite propagator.
  - `gs::GS`: Ground station coordinates `(lat, lon [, height])` [rad, m].
  - `sat1::OrbitPropagatorSgp4`: First satellite propagator.
  - `sat2::OrbitPropagatorSgp4`: Second satellite propagator.

# Keyword Arguments

  - `transmitter_diameter_m::Real=0.6`: Transmitting telescope diameter [m].
  - `receiver_diameter_m::Real=0.6`: Receiving telescope diameter [m].
  - `wavelength_nm::Real=1550`: Wavelength [nm].
  - `conditions::Conditions=clear`: Atmospheric conditions.
  - `time::Time=missing`: Evaluation time in Julian days (Real) or DateTime
    (default: satellite epoch).
  - `min_angle::Real=deg2rad(20)`: Minimum elevation angle [rad] (sat-gs only).

# Returns

  - `Union{FreespaceChannel, Nothing}`: `FreespaceChannel` instance or `nothing` if
    visibility constraint not met.

# Notes

  - See [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)
    for reference of `OrbitPropagatorSgp4`.
"""
struct FreespaceChannel <: AbstractChannel
    distance_m::Float64
    elevation_angle_rad::Float64
    min_altitude_m::Float64
    transmitter_diameter_m::Float64
    receiver_diameter_m::Float64
    wavelength_nm::Float64
    conditions::Conditions
    light_condition::LightCondition

    FreespaceChannel(
        distance_m::Real,
        elevation_angle_rad::Real;
        min_altitude_m::Real=0,
        transmitter_diameter_m::Real=0.6,
        receiver_diameter_m::Real=0.6,
        wavelength_nm::Real=1550,
        conditions::Conditions=clear,
        light_condition::LightCondition=umbra,
    ) = new(
        Float64(distance_m),
        Float64(elevation_angle_rad),
        Float64(min_altitude_m),
        Float64(transmitter_diameter_m),
        Float64(receiver_diameter_m),
        Float64(wavelength_nm),
        conditions,
        light_condition,
    )
end
function FreespaceChannel(
    sat::OrbitPropagatorSgp4,
    gs::GS;
    transmitter_diameter_m::Real=0.6,
    receiver_diameter_m::Real=0.6,
    wavelength_nm::Real=1550,
    time::Time=missing,
    conditions::Conditions=clear,
    min_angle::Real=deg2rad(20),
)
    gs_height = length(gs) > 2 ? gs[3] : 0

    time = _get_time_JD(time, sat)
    Δt_s = _JD_to_seconds(time - sat.sgp4d.epoch)

    sat_sv = Propagators.propagate!(sat, Δt_s, OrbitStateVector)
    sat_pos = eci_to_ecef(sat_sv.r; time=time)

    if !is_ground_facility_visible(sat_pos, gs[1], gs[2], gs_height, min_angle)
        return nothing
    end

    n, e, d = ecef_to_ned(sat_pos, gs[1], gs[2], gs_height; translate=true)
    elevation_angle_rad = atan(-d / √(n^2 + e^2))

    min_altitude_m = d < 0 ? gs_height : norm(sat_pos)
    distance_m = norm(geodetic_to_ecef(gs) .- sat_pos)

    light_condition = light_at_point(sat_sv.r, time)

    return FreespaceChannel(
        distance_m,
        elevation_angle_rad;
        min_altitude_m,
        transmitter_diameter_m,
        receiver_diameter_m,
        wavelength_nm,
        conditions,
        light_condition,
    )
end
function FreespaceChannel(
    sat1::OrbitPropagatorSgp4,
    sat2::OrbitPropagatorSgp4;
    transmitter_diameter_m::Real=0.6,
    receiver_diameter_m::Real=0.6,
    wavelength_nm::Real=1550,
    time::Time=missing,
    conditions::Conditions=clear,
)
    time = _get_time_JD(time, sat1, sat2)

    Δt1_s = _JD_to_seconds(time - sat1.sgp4d.epoch)
    sat1_sv = Propagators.propagate!(sat1, Δt1_s, OrbitStateVector)
    sat1_pos = eci_to_ecef(sat1_sv.r; time=time)

    Δt2_s = _JD_to_seconds(time - sat2.sgp4d.epoch)
    sat2_sv = Propagators.propagate!(sat2, Δt2_s, OrbitStateVector)
    sat2_pos = eci_to_ecef(sat2_sv.r; time=time)

    intersection_points = ellipsoid_line_intersection(
        SEMIMAJOR_RADIUS, SEMIMINOR_RADIUS, sat1_pos, sat2_pos; only_between=true
    )

    if length(intersection_points) > 0
        return nothing
    end

    min_altitude_m = min(norm(sat1_pos), norm(sat2_pos))
    distance_m = norm(sat1_pos .- sat2_pos)

    light_condition = min(light_at_point(sat1_sv.r, time), light_at_point(sat2_sv.r, time))

    return FreespaceChannel(
        distance_m,
        0;
        min_altitude_m,
        transmitter_diameter_m,
        receiver_diameter_m,
        wavelength_nm,
        conditions,
        light_condition,
    )
end
