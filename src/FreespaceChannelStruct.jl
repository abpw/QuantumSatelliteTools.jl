using ..AstronomyGeometry
using Dates: DateTime, datetime2julian
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using SatelliteToolboxBase: OrbitStateVector
using SatelliteAnalysis: is_ground_facility_visible
using SatelliteToolboxTransformations: ecef_to_ned

@enum Conditions begin
    clear
    fog
    rain
    snow
end
# TODO doctring for struct

struct FreespaceChannel
    distance_m::Number
    elevation_angle_rad::Number
    min_altitude_m::Number
    transmitter_diameter_m::Number
    receiver_diameter_m::Number
    wavelength_nm::Number
    conditions::Conditions
    light_condition::LightCondition
end

"""
Create a FreespaceChannel with specified parameters.

# Arguments
- `distance_m::Number`: Distance in meters.
- `elevation_angle_rad::Number`: Elevation angle in radians.

# Keyword Arguments
- `transmitter_diameter_m::Number`: Diameter of the transmitting telescope in meters (default 0.6).
- `receiver_diameter_m::Number`: Diameter of the receiving telescope in meters (default 0.6).
- `wavelength_nm::Number`: Wavelength in nanometers (default 1550).
- `min_altitude_m::Number`: Minimum altitude in meters (default 0).
- `conditions::Conditions`: Atmospheric conditions (default clear).
- `light_condition::LightCondition`: Lighting condition (default umbra).

# Returns
- `FreespaceChannel`: A new FreespaceChannel instance.
"""
function FreespaceChannel(distance_m::Number, elevation_angle_rad::Number; min_altitude_m::Number=0, transmitter_diameter_m::Number=0.6, receiver_diameter_m::Number=0.6, wavelength_nm::Number=1550, conditions::Conditions=clear, light_condition::LightCondition=umbra)
    return FreespaceChannel(distance_m, elevation_angle_rad, min_altitude_m, transmitter_diameter_m, receiver_diameter_m, wavelength_nm, conditions, light_condition)
end

"""
Create a FreespaceChannel between a satellite and a ground station.

# Arguments
- `sat::OrbitPropagatorSgp4{Float64, Float64}`: Satellite propagator (see [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)).
- `gs::GS`: Ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).

# Keyword Arguments
- `transmitter_diameter_m::Number`: Diameter of the transmitting telescope in meters (default 0.6).
- `receiver_diameter_m::Number`: Diameter of the receiving telescope in meters (default 0.6).
- `wavelength_nm::Number`: Wavelength in nanometers (default 1550).
- `time::Union{Number, DateTime}`: Julian time of calculation (defaults to satellite epoch).
- `conditions::Conditions`: Atmospheric conditions (default clear).
- `min_θ`: Minimum elevation angle in radians (default 20 degrees).

# Returns
- `FreespaceChannel` instance or `nothing` if ground station is not visible.
"""
function FreespaceChannel(sat::OrbitPropagatorSgp4{Float64, Float64}, gs::GS; transmitter_diameter_m::Number=60, receiver_diameter_m::Number=60, wavelength_nm::Number=1550, time::Union{Number,DateTime,Missing}=missing, conditions::Conditions=clear, min_θ=deg2rad(20))
    if time === missing
        time = sat.sgp4d.epoch
    end
    if typeof(time) == DateTime
        time = datetime2julian(time)
    end
    gs_height = length(gs) > 2 ? gs[3] : 0
    sat_sv = Propagators.propagate!(sat, time - sat.sgp4d.epoch, OrbitStateVector)
    sat_pos = eci_to_ecef(sat_sv.r, time=time)
    if !is_ground_facility_visible(sat_pos, gs[1], gs[2], gs_height, min_θ)
        return nothing
    end
    gs_height = 0
    if length(gs) > 2
        gs_height = gs[3]
    end
    n, e, d = ecef_to_ned(sat_pos, gs[1], gs[2], gs_height, translate=true)
    elevation_angle_rad = atan(-d / √(n^2 + e^2))
    min_altitude_m = d < 0 ? gs_height : √(sat_pos[1]^2 + sat_pos[2]^2 + sat_pos[3]^2)
    light_condition = light_at_point(sat_sv.r, time)
    distance_m = norm(geodetic_to_ecef(gs) .- sat_pos)
    return FreespaceChannel(distance_m, elevation_angle_rad, min_altitude_m, transmitter_diameter_m, receiver_diameter_m, wavelength_nm, conditions, light_condition)
end

"""
Create a FreespaceChannel between two satellites.

# Arguments
- `sat1::OrbitPropagatorSgp4{Float64, Float64}`: Satellite propagator (see [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)).
- `sat2::OrbitPropagatorSgp4{Float64, Float64}`: Satellite propagator.

# Keyword Arguments
- `transmitter_diameter_m::Number`: Diameter of the transmitting telescope in meters (default 0.6).
- `receiver_diameter_m::Number`: Diameter of the receiving telescope in meters (default 0.6).
- `wavelength_nm::Number`: Wavelength in nanometers (default 1550).
- `time::Union{Number, DateTime}`: Julian time of calculation (defaults to satellite epoch).
- `conditions::Conditions`: Atmospheric conditions (default clear).

# Returns
- `FreespaceChannel` instance or `nothing` if satellites are not visible to each other.
"""
function FreespaceChannel(sat1::OrbitPropagatorSgp4{Float64, Float64}, sat2::OrbitPropagatorSgp4{Float64, Float64}; transmitter_diameter_m::Number=60, receiver_diameter_m::Number=60, wavelength_nm::Number=1550, time::Union{Number,DateTime,Missing}=missing, conditions::Conditions=clear)
    if time === missing
        time = sat.sgp4d.epoch
    end
    if typeof(time) == DateTime
        time = datetime2julian(time)
    end

    sat1_sv = Propagators.propagate!(sat1, time - sat1.sgp4d.epoch, OrbitStateVector)
    sat1_pos = eci_to_ecef(sat1_sv.r, time=time)

    sat2_sv = Propagators.propagate!(sat2, time - sat2.sgp4d.epoch, OrbitStateVector)
    sat2_pos = eci_to_ecef(sat2_sv.r, time=time)

    visibility_check = ellipsoid_line_intersection(semimajor_radius, semiminor_radius, sat1_pos, sat2_pos, only_between=true)

    if length(visibility_check) > 0
        return nothing
    end

    min_altitude_m = min(√(sat1_pos[1]^2 + sat1_pos[2]^2 + sat1_pos[3]^2), √(sat2_pos[1]^2 + sat2_pos[2]^2 + sat2_pos[3]^2))
    # TODO light at worst case of both points
    light_condition = light_at_point(sat1_sv.r, time)

    distance_m = norm(sat1_pos .- sat2_pos)

    FreespaceChannel(distance_m, 0, min_altitude_m, transmitter_diameter_m, receiver_diameter_m, wavelength_nm, conditions, light_condition)
end
