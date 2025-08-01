using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using SatelliteToolboxTle: TLE
using SatelliteToolboxTransformations
using SatelliteToolboxTransformations: WGS84_ELLIPSOID, Ellipsoid
import SatelliteToolboxTransformations: ecef_to_geodetic, geodetic_to_ecef
using Dates: DateTime, datetime2julian
using SatelliteAnalysis: is_ground_facility_visible, lighting_condition
using SatelliteToolboxBase: JD_J2000
using StaticArrays: SVector

# of Earth, in meters
const semimajor_radius = 6378137.0
const semiminor_radius = 6356752.3

const G = 6.6743e-11
const earth_mass_kg = 5.9722e24
const seconds_per_day = 60*60*24

const equatorial_circumference_km = semimajor_radius / 500 * π
const sin_60 = √3/2

GS = Union{NTuple{2, Number}, NTuple{3, Number}, SVector{2, Number}, SVector{3, Number}}
Num64 = Union{Float64, Int64}
Point3D = Union{SVector{3, Float64}, SVector{3, Int64}, NTuple{3, Num64}}

"""
Compute the Euclidean norm of a 3D point.

# Arguments
- `point::Point3D`: A 3D point given as a tuple or static vector.

# Returns
- `Float64`: The Euclidean distance from the origin.
"""
function norm(point::Point3D)
    √sum(point.^2)
end

"""
Convert geodetic coordinates to ECEF.

# Arguments
- `gs::GS`: Ground station coordinates (geodetic latitude in radians, geodetic longitude in radians[, height in meters]).

# Keyword Arguments
- `ellipsoid`: Ellipsoid model to use (defaults to WGS84, see [SatelliteToolboxTransformations.jl](https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl)).

# Returns
- ECEF coordinates in meters.
"""
function geodetic_to_ecef(gs::GS; ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID) where T<:Number
    gs_height = 0
    if length(gs) > 2
        gs_height = gs[3]
    end
    return geodetic_to_ecef(gs[1], gs[2], gs_height, ellipsoid=ellipsoid)
end

"""
Convert ECEF coordinates (x, y, z) to geodetic coordinates.

# Arguments
- `x::Num64`: X coordinate in ECEF.
- `y::Num64`: Y coordinate in ECEF.
- `z::Num64`: Z coordinate in ECEF.

# Keyword Arguments
- `ellipsoid`: Ellipsoid model to use (defaults to WGS84).

# Returns
- Geodetic coordinates (geodetic latitude in radians, geodetic latitude in radians, height in meters).
"""
function ecef_to_geodetic(x::Num64, y::Num64, z::Num64; ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID) where T<:Number
    return ecef_to_geodetic([x, y, z], ellipsoid=ellipsoid)
end

"""
Convert ECEF coordinates given as a 3D point to geodetic coordinates.

# Arguments
- `point_ecef::Point3D`: ECEF coordinates in meters.

# Keyword Arguments
- `ellipsoid`: Ellipsoid model to use (defaults to WGS84).

# Returns
- Geodetic coordinates (geodetic latitude in radians, geodetic latitude in radians, height in meters).
"""
function ecef_to_geodetic(point_ecef::Point3D; ellipsoid::Ellipsoid{T} = WGS84_ELLIPSOID) where T<:Number
    return ecef_to_geodetic([point_ecef[1], point_ecef[2], point_ecef[3]], ellipsoid=ellipsoid)
end

"""
Convert ECI coordinates to ECEF coordinates.

# Arguments
- `sat_sv::OrbitStateVector`: Satellite state vector (see [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl)).

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of conversion (defaults to state vector time).

# Returns
- ECEF position vector in meters.
"""
function eci_to_ecef(sat_sv::OrbitStateVector; time::Union{Number, DateTime, Missing}=missing)
    if time === missing
        time = sat_sv.t
    end
    if typeof(time) == DateTime
        time = datetime2julian(time)
    end
    return SatelliteToolboxTransformations.sv_eci_to_ecef(sat_sv, SatelliteToolboxTransformations.TEME(), SatelliteToolboxTransformations.PEF(), time).r
end

"""
Convert ECI satellite propagator to ECEF coordinates.

# Arguments
- `sat::OrbitPropagatorSgp4`: Satellite propagator.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Time of conversion (defaults to satellite epoch).

# Returns
- ECEF position vector.
"""
function eci_to_ecef(sat::OrbitPropagatorSgp4; time::Union{Number, DateTime, Missing}=missing)
    return eci_to_ecef(Propagators.propagate!(sat, 0), time=time)
end

"""
Convert a 3D ECI point to ECEF coordinates.

# Arguments
- `point_eci::Point3D`: ECI 3D point in meters.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of conversion.

# Returns
- ECEF position vector in meters.
"""
function eci_to_ecef(point_eci::Point3D; time::Union{Number, DateTime})
    return eci_to_ecef(OrbitStateVector(0, [point_eci[1], point_eci[2], point_eci[3]], [0,0,0]), time=time)
end

"""
Convert ECI coordinates (x, y, z) to ECEF coordinates.

# Arguments
- `x::Num64`: ECI X coordinate in meters.
- `y::Num64`: ECI Y coordinate in meters.
- `z::Num64`: ECI Z coordinate in meters.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of conversion.

# Returns
- ECEF position vector in meters.
"""
function eci_to_ecef(x::Num64, y::Num64, z::Num64; time::Union{Number, DateTime})
    return eci_to_ecef((x, y, z), time=time)
end

"""
Convert ECI satellite state vector to geodetic coordinates.

# Arguments
- `sat_sv::OrbitStateVector`: Satellite state vector (see [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl)).

# Keyword Arguments
- `time::Union{Number, DateTime}`: Time of conversion (defaults to state vector time).

# Returns
- Geodetic coordinates (geodetic latitude in radians, geodetic lon in radians, height in meters).
"""
function eci_to_geodetic(sat_sv::OrbitStateVector; time::Union{Number, DateTime, Missing}=missing)
    if time === missing
        time = sat_sv.t
    end
    return ecef_to_geodetic(eci_to_ecef(sat_sv, time=time))
end

"""
Convert satellite propagator to geodetic coordinates.

# Arguments
- `sat::OrbitPropagatorSgp4`: Satellite propagator (see [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)).

# Keyword Arguments
- `time::Union{Number, DateTime}`: Time of conversion (defaults to the satellite epoch).

# Returns
- Geodetic coordinates (geodetic latitude in radians, geodetic lon in radians, height in meters).
"""
function eci_to_geodetic(sat::OrbitPropagatorSgp4; time::Union{Number, DateTime, Missing}=missing)
    return eci_to_geodetic(Propagators.propagate!(sat, 0, OrbitStateVector), time=time)
end

"""
Convert a 3D ECI point to geodetic coordinates.

# Arguments
- `point_eci::Point3D`: 3D point in ECI coordinates.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of conversion.

# Returns
- Geodetic coordinates (geodetic latitude in radians, geodetic lon in radians, height in meters).
"""
function eci_to_geodetic(point_eci::Point3D; time::Union{Number, DateTime})
    return eci_to_geodetic(OrbitStateVector(0, [point_eci[1], point_eci[2], point_eci[3]], [0,0,0]), time=time)
end

"""
Convert ECI coordinates (x, y, z) to geodetic coordinates.

# Arguments
- `x::Num64`: ECI X coordinate in meters.
- `y::Num64`: ECI Y coordinate in meters.
- `z::Num64`: ECI Z coordinate in meters.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of conversion.

# Returns
- Geodetic coordinates (latitude in radians, geodetic lon in radians, height in meters).
"""
function eci_to_geodetic(x::Num64, y::Num64, z::Num64; time::Union{Number, DateTime})
    return eci_to_geodetic((x, y, z), time=time)
end

"""
Compute the intersection points between an ellipsoid and a line segment.

# Arguments
- `major::Num64`: Semimajor axis of the ellipsoid.
- `minor::Num64`: Semiminor axis of the ellipsoid.
- `point1::Point3D`: First point defining the line in the ECEF or ECI reference frame.
- `point2::Point3D`: Second point defining the line in the ECEF or ECI reference frame.

# Keyword Arguments
- `only_between::Bool`: If true, only return intersections between point1 and point2.

# Returns
- A vector of 0, 1, or 2 intersection points as tuples in the same reference frame as the inputs (ECI or ECEF).
"""
function ellipsoid_line_intersection(major::Num64, minor::Num64, point1::Point3D, point2::Point3D; only_between::Bool=false)
    x1, y1, z1 = point1
    x2, y2, z2 = point2

    # coefficients for quadratic equation describing the intersection
    # of the line between the ground station and the satellite and
    # the atmosphere ellipsoid
    c = (x1/major)^2 + (y1/major)^2 + (z1/minor)^2 - 1
    b = 2*x1*(x2-x1)/major^2
    b += 2*y1*(y2-y1)/major^2
    b += 2*z1*(z2-z1)/minor^2
    a = ((x2-x1)/major)^2 + ((y2-y1)/major)^2 + ((z2-z1)/minor)^2

    # solving the quadratic equation gives the intersection parameters
    # for the line between the ground station and the satellite

    # no real solutions case (no intersection)
    if b^2-4*a*c < 0
        return []
    end
    
    out = []

    # two solutions (intersects twice)
    if b^2-4*a*c > 0
        t2 = (-1*b - sqrt(b^2-4*a*c))/2/a
        out_x2 = x1 + t2*(x2 - x1); out_y2 = y1 + t2*(y2 - y1); out_z2 = z1 + t2*(z2 - z1)
        out2 = (out_x2, out_y2, out_z2)
        if !only_between || (norm(point1 .- out2) + norm(out2 .- point2) ≈ norm(point1 .- point2) && !(norm(point1 .- out2) ≈ 0) && !(norm(point2 .- out2) ≈ 0))
            push!(out, out2)
        end
    end

    # one or two solutions
    t1 = (-1*b + sqrt(b^2-4*a*c))/2/a
    out_x1 = x1 + t1*(x2 - x1); out_y1 = y1 + t1*(y2 - y1); out_z1 = z1 + t1*(z2 - z1)
    out1 = (out_x1, out_y1, out_z1)
    if !only_between || (norm(point1 .- out1) + norm(out1 .- point2) ≈ norm(point1 .- point2) && !(norm(point1 .- out1) ≈ 0) && !(norm(point2 .- out1) ≈ 0))
        push!(out, out1)
    end

    return out
end

"""
Compute the distance between two ground stations given in geodetic coordinates.

# Arguments
- `gs1::GS`: First ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).
- `gs2::GS`: Second ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).

# Returns
- Distance between the ground stations in meters. Uses Euclidean distance if there is line of sight (one ground station must have positive height). Otherwise uses the [haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).
"""
function gs_gs_distance(gs1::GS, gs2::GS)
    if length(ellipsoid_line_intersection(semimajor_radius, semiminor_radius, geodetic_to_ecef(gs1), geodetic_to_ecef(gs2), only_between=true)) == 0
        # use line of sight distance if there's line of sight between them
        return norm(geodetic_to_ecef(gs1) .- geodetic_to_ecef(gs2))
    end
    #height1 = 0 ? length(gs1) < 3 : gs1[3]
    #height2 = 0 ? length(gs2) < 3 : gs2[3]
    lat1, lon1 = gs1[1:2]; lat2, lon2 = gs2[1:2]
    hav = 1/2 - cos(lat2 - lat1)/2 + cos(lat1)*cos(lat2)*(1/2 - cos(lon2 - lon1)/2)
    theta = acos(1-hav*2)
    return theta * (semimajor_radius + semiminor_radius)/2
end

"""
Compute the distance between two ground stations given in ECEF coordinates.

# Arguments
- `::Val{:ECEF}`: Val symbol to specify ECEF input.
- `gs1::GS`: First ground station ECEF coordinates in meters.
- `gs2::GS`: Second ground station ECEF coordinates in meters.

# Returns
- Distance in meters.
"""
function gs_gs_distance(::Val{:ECEF}, gs1::GS, gs2::GS)
    return gs_gs_distance(ecef_to_geodetic(gs1), ecef_to_geodetic(gs2))
end

"""
Compute the distance between two ground stations given in geodetic coordinates.

# Arguments
- `::Val{:lat_lon}`: Val symbol to specify lat/lon input.
- `gs1::GS`: First ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).
- `gs2::GS`: Second ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).

# Returns
- Distance in meters.
"""
function gs_gs_distance(::Val{:lat_lon}, gs1::GS, gs2::GS)
    return gs_gs_distance(gs1, gs2)
end

"""
Compute the distance between two satellites at a given time.

# Arguments
- `sat_prop1::OrbitPropagatorSgp4`: First satellite propagator (see [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)).
- `sat_prop2::OrbitPropagatorSgp4`: Second satellite propagator.

# Keyword Arguments
- `time::Union{Number, DateTime}`: Time of distance calculation (defaults to the later epoch of the satellites).

# Returns
- Distance in meters or `nothing` if obstructed.
"""
function sat_sat_distance(sat_prop1::OrbitPropagatorSgp4, sat_prop2::OrbitPropagatorSgp4; time::Union{Number, DateTime, Missing}=missing)
    if time === missing
        time = max(sat_prop1.sgp4d.epoch, sat_prop2.sgp4d.epoch)
    end
    pos1 = Tuple(Propagators.propagate!(sat_prop1, time - sat_prop1.sgp4d.epoch)[1])
    pos2 = Tuple(Propagators.propagate!(sat_prop2, time - sat_prop1.sgp4d.epoch)[1])
    visibility_check = ellipsoid_line_intersection(semimajor_radius, semiminor_radius, pos1, pos2)
    if length(visibility_check) > 0
        return nothing
    end
    return √((pos1[1]-pos2[1])^2 + (pos1[2]-pos2[2])^2 + (pos1[3]-pos2[3])^2)
end

"""
Compute the distance between two satellites given their TLEs.

# Arguments
- `sat_tle1::TLE`: First satellite TLE (see GenerateTLEs or [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)).
- `sat_tle2::TLE`: Second satellite TLE (see GenerateTLEs or [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)).

# Keyword Arguments
- `time::Union{Number, DateTime}`: Time of distance calculation (defaults to the later epoch of the satellites).

# Returns
- Distance in meters or `nothing` if obstructed.
"""
function sat_sat_distance(sat_tle1::TLE, sat_tle2::TLE; time::Union{Number, DateTime, Missing}=missing)
    return sat_sat_distance(Propagators.init(Val(:SGP4), sat_tle1), Propagators.init(Val(:SGP4), sat_tle2), time=time)    
end

@enum LightCondition begin
    sunlight=true
    penumbra
    umbra=false
end

"""
Calculate the position of the sun in ECI coordinates at a given Julian time.

# Arguments
- `time::Float64`: Julian time of evaluation.

# Returns
- ECI Sun position vector in meters.
"""
function sun_position(time::Float64)
    d = time - JD_J2000
    # https://astronomy.stackexchange.com/questions/28802/calculating-the-sun-s-position-in-eci
    # Calculate parameters
    L = (280.4606184 + ((36000.77005361 / 36525) * d))*π/180 # mean longitude, in radians
    g = (357.5277233 + ((35999.05034 / 36525) * d))*π/180 # mean anomaly, in radians
    p = L + ((1.914666471 * sin(g)) + (0.918994643 * sin(2*g)))*π/180 # ecliptic longitude lambda, in radians
    q = (23.43929 - ((46.8093/3600) * (d / 36525)))*π/180 # obliquity of ecliptic plane ϵ, in radians

    # Calculate unit directional vector in ECI coordinates
    direction_vector = [cos(p), cos(q) * sin(p), sin(q) * sin(p)]

    # Calculate distance to sun and scale the unit vector
    a = 1.000140612 - (0.016708617 * cos(g)) - (0.000139589 * cos(2*g)) # distance from Earth's center to Sun's center in astronomical units (AU)
    m = a * 149597870700 # center-to-center distance from Earth to Sun in meters
    return m .* direction_vector # distance to sun in meters
end

"""
Calculate the position of the sun in ECI coordinates at a given DateTime.

# Arguments
- `time::DateTime`: Julian time of evaluation.

# Returns
- ECI Sun position vector in meters.
"""
function sun_position(time::DateTime)
    return sun_position(datetime2julian(time))
end

"""
Determine the lighting condition at a given ECI point and time.

# Arguments
- `point_eci::Point3D`: ECI point in meters.
- `time::Union{Number, DateTime}`: Julian time of evaluation.

# Returns
- `LightCondition`: The lighting condition (sunlight, penumbra, umbra).
"""
function light_at_point(point_eci::Point3D, time::Union{Number, DateTime})
    sun_pos = sun_position(time)
    light_symbol = lighting_condition([point_eci[1], point_eci[2], point_eci[3]], sun_pos)
    if light_symbol == :sunlight
        return sunlight
    elseif light_symbol == :penumbra
        return penumbra
    elseif light_symbol == :umbra
        return umbra
    end
end

@enum Conditions begin 
    clear
    fog
    rain
    snow
end

struct FreespaceChannel
    distance_m::Number
    elevation_angle_rad::Number
    min_altitude_m::Number
    conditions::Conditions
    light_condition::LightCondition
end

"""
Create a FreespaceChannel with specified parameters.

# Arguments
- `distance_m::Num64`: Distance in meters.
- `elevation_angle_rad::Num64`: Elevation angle in radians.

# Keyword Arguments
- `min_altitude_m::Num64`: Minimum altitude in meters (default 0).
- `conditions::Conditions`: Atmospheric conditions (default clear).
- `light_condition::LightCondition`: Lighting condition (default umbra).

# Returns
- `FreespaceChannel`: A new FreespaceChannel instance.
"""
function FreespaceChannel(distance_m::Num64, elevation_angle_rad::Num64; min_altitude_m::Num64=0, conditions::Conditions=clear, light_condition::LightCondition=umbra)
    return FreespaceChannel(distance_m, elevation_angle_rad, min_altitude_m, conditions, light_condition)    
end
function FreespaceChannel(distance_m::Num64, elevation_angle_rad::Num64; min_altitude_m::Num64=0, conditions::Conditions=clear, light_condition::LightCondition=umbra)
    return FreespaceChannel(distance_m, elevation_angle_rad, min_altitude_m, conditions, light_condition)    
end

"""
Create a FreespaceChannel between a satellite and a ground station.

# Arguments
- `sat::OrbitPropagatorSgp4`: Satellite propagator (see [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)).
- `gs::GS`: Ground station coordinates (geodetic latitude in radians, geodetic lon in radians[, height in meters]).

# Keyword Arguments
- `time::Union{Number, DateTime}`: Julian time of calculation (defaults to satellite epoch).
- `conditions::Conditions`: Atmospheric conditions (default clear).
- `min_θ`: Minimum elevation angle in radians (default 20 degrees).

# Returns
- `FreespaceChannel` instance or `nothing` if ground station is not visible.
"""
function FreespaceChannel(sat::OrbitPropagatorSgp4, gs::GS; time::Union{Number, DateTime, Missing}=missing, conditions::Conditions=clear, min_θ=deg2rad(20))
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
    n, e, d = SatelliteToolboxTransformations.ecef_to_ned(sat_pos, gs[1], gs[2], gs_height, translate=true)
    elevation_angle_rad = atan(-d/√(n^2+e^2))
    min_altitude_m = d < 0 ? gs_height : √(sat_pos[1]^2 + sat_pos[2]^2 + sat_pos[3]^2)
    light_condition = light_at_point(sat_sv.r, time)
    distance_m = norm(geodetic_to_ecef(gs).-sat_pos)
    return FreespaceChannel(distance_m, elevation_angle_rad, min_altitude_m, conditions, light_condition)
end