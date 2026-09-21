using Dates: DateTime, datetime2julian

import SatelliteToolboxTransformations: ecef_to_geodetic, geodetic_to_ecef

using SatelliteAnalysis: lighting_condition
using SatelliteToolboxBase: JD_J2000, WGS84_ELLIPSOID, Ellipsoid, OrbitStateVector
using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
using SatelliteToolboxTle: TLE
using SatelliteToolboxTransformations: PEF, TEME, sv_eci_to_ecef

using ..constants: SEMIMAJOR_RADIUS, SEMIMINOR_RADIUS
using ..helpers: norm, _to_float64, validate_coordinates, _JD_to_seconds, _get_time_JD
using ..types: Point3D, InputGS, GS, Time, AbstractChannel, sunlight, penumbra, umbra

"""
    geodetic_to_ecef(
        gs::InputGS;
        ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID
    ) where {T<:Real} -> Point3D

Convert geodetic coordinates `gs` [rad, m], to the Earth-Centered
Earth-Fixed (ECEF) reference frame with respect to the specified `ellipsoid` model.

See [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl/)
for reference of WGS84_ELLIPSOID.

# Returns

  - `Point3D`: ECEF coordinates [m].

# Throws

  - `ArgumentError`: If any of the ground station coordinates are not finite numbers.

# See Also

  - [`ecef_to_geodetic`](#ecef_to_geodetic)
"""
function geodetic_to_ecef(gs::InputGS; ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID) where {T<:Real}
    validate_coordinates(gs)

    gs_height = length(gs) > 2 ? gs[3] : 0

    return geodetic_to_ecef(gs[1], gs[2], gs_height; ellipsoid=ellipsoid)
end

"""
    ecef_to_geodetic(
        x::Real, y::Real, z::Real;
        ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID
    ) where {T<:Real} -> Point3D
    ecef_to_geodetic(
        point_ecef::Point3D;
        ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID
    ) where {T<:Real} -> Point3D

Convert Earth-Centered Earth-Fixed (ECEF) coordinates, given by `x, y, z` [m] or
`point_ecef` [m], to the geodetic reference frame with respect to the specified
`ellipsoid` model.

See [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl/)
for reference of WGS84_ELLIPSOID.

# Returns

  - `Point3D`: Geodetic coordinates `(latitude, longitude, height)` [rad, m].

# Throws

  - `ArgumentError`: If any of the ECEF coordinates are not finite numbers.

# See Also

  - [`geodetic_to_ecef`](#geodetic_to_ecef)
"""
function ecef_to_geodetic(
    x::Real, y::Real, z::Real; ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID
) where {T<:Real}
    validate_coordinates((x, y, z))

    return ecef_to_geodetic([x, y, z]; ellipsoid=ellipsoid)
end
function ecef_to_geodetic(
    point_ecef::Point3D; ellipsoid::Ellipsoid{T}=WGS84_ELLIPSOID
) where {T<:Real}
    validate_coordinates(point_ecef)

    return ecef_to_geodetic(
        [point_ecef[1], point_ecef[2], point_ecef[3]]; ellipsoid=ellipsoid
    )
end

"""
    eci_to_ecef(sat_sv::OrbitStateVector; time::Time=missing) -> Point3D
    eci_to_ecef(sat::OrbitPropagatorSgp4; time::Time=missing) -> Point3D
    eci_to_ecef(point_eci::Point3D; time::Time=missing) -> Point3D
    eci_to_ecef(x::Real, y::Real, z::Real; time::Time=missing) -> Point3D

Convert Earth-Centered Inertial (ECI) coordinates, given by `sat_sv`, `sat`,
`point_eci` [m], or `x, y, z` [m], to the Earth-Centered Earth-Fixed (ECEF) reference frame.

If `time` [Julian days or DateTime] is not provided, the epoch of the state vector or
propagator will be used.

See [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl) for
reference of `OrbitStateVector` and
[SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)
for reference of `OrbitPropagatorSgp4`.

# Returns

  - `Point3D`: ECEF coordinates [m].

# Throws

  - `ArgumentError`: If any of the ECI coordinates are not finite numbers.
"""
function eci_to_ecef(sat_sv::OrbitStateVector; time::Time=missing)
    validate_coordinates(sat_sv.r)

    time = _get_time_JD(time, sat_sv)

    return sv_eci_to_ecef(sat_sv, TEME(), PEF(), time).r
end
function eci_to_ecef(sat::OrbitPropagatorSgp4; time::Time=missing)
    return eci_to_ecef(Propagators.propagate!(sat, 0, OrbitStateVector); time=time)
end
function eci_to_ecef(point_eci::Point3D; time::Time=missing)
    return eci_to_ecef(
        OrbitStateVector(0, [point_eci[1], point_eci[2], point_eci[3]], [0, 0, 0]);
        time=time,
    )
end
function eci_to_ecef(x::Real, y::Real, z::Real; time::Time=missing)
    return eci_to_ecef((x, y, z); time=time)
end

"""
    eci_to_geodetic(sat_sv::OrbitStateVector; time::Time=missing) -> Point3D
    eci_to_geodetic(sat::OrbitPropagatorSgp4; time::Time=missing) -> Point3D
    eci_to_geodetic(point_eci::Point3D; time::Time=missing) -> Point3D
    eci_to_geodetic(x::Real, y::Real, z::Real; time::Time=missing) -> Point3D

Convert Earth-Centered Inertial (ECI) coordinates, given by `sat_sv`, `sat`,
`point_eci` [m], or `x, y, z` [m], to the geodetic reference frame.

If `time` [Julian days or DateTime] is not provided, the epoch of the state vector or
propagator will be used.

See [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl) for
reference of `OrbitStateVector` and
[SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)
for reference of `OrbitPropagatorSgp4`.

# Returns

  - `Point3D`: Geodetic coordinates `(latitude, longitude, height)` [rad, m].
"""
function eci_to_geodetic(sat_sv::OrbitStateVector; time::Time=missing)
    validate_coordinates(sat_sv.r)

    time = _get_time_JD(time, sat_sv)

    return ecef_to_geodetic(eci_to_ecef(sat_sv; time=time))
end
function eci_to_geodetic(sat::OrbitPropagatorSgp4; time::Time=missing)
    return eci_to_geodetic(Propagators.propagate!(sat, 0, OrbitStateVector); time=time)
end
function eci_to_geodetic(point_eci::Point3D; time::Time=missing)
    return eci_to_geodetic(
        OrbitStateVector(0, [point_eci[1], point_eci[2], point_eci[3]], [0, 0, 0]);
        time=time,
    )
end
function eci_to_geodetic(x::Real, y::Real, z::Real; time::Time=missing)
    return eci_to_geodetic((x, y, z); time=time)
end

"""
    ellipsoid_line_intersection(
        major::Real, minor::Real, point1::Point3D, point2::Point3D;
        only_between::Bool=false
    ) -> Vector{Point3D}

Compute intersection points between an ellipsoid and a line segment.

# Arguments

  - `major::Real`: Semimajor axis of ellipsoid [m].
  - `minor::Real`: Semiminor axis of ellipsoid [m].
  - `point1::Point3D`: First point of line segment [m].
  - `point2::Point3D`: Second point of line segment [m].

# Keyword Arguments

  - `only_between::Bool=false`: If `true`, only return intersection points that lie
    between `point1` and `point2`.

# Returns

  - `Vector{Point3D}`: 0 - 2 intersection points
    (same reference frame as inputs; ECI or ECEF).

# Notes
  - `point1` and `point2` must be in the same reference frame (ECI or ECEF).
"""
function ellipsoid_line_intersection(
    major::Real, minor::Real, point1::Point3D, point2::Point3D; only_between::Bool=false
)
    if !isfinite(major) || !isfinite(minor)
        throw(ArgumentError("Ellipsoid axes must be finite, got $(major), $(minor)."))
    end
    validate_coordinates(point1)
    validate_coordinates(point2)

    x₁, y₁, z₁ = point1
    x₂, y₂, z₂ = point2

    # quadratic coefficients for line-ellipsoid intersection
    c = (x₁ / major)^2 + (y₁ / major)^2 + (z₁ / minor)^2 - 1
    b =
        2 * x₁ * (x₂ - x₁) / major^2 +
            2 * y₁ * (y₂ - y₁) / major^2 +
            2 * z₁ * (z₂ - z₁) / minor^2
    a = ((x₂ - x₁) / major)^2 + ((y₂ - y₁) / major)^2 + ((z₂ - z₁) / minor)^2

    Δ = b^2 - 4 * a * c

    # no real solutions case (no intersection)
    if Δ < 0
        return []
    end

    # one or two solutions case (intersects once or twice)
    t₁ = (-b - √Δ) / 2a
    t₂ = (-b + √Δ) / 2a

    t_values = isapprox(t₁, t₂) ? [t₁] : sort([t₁, t₂])

    out = []
    for t in t_values
        point = (x₁ + t * (x₂ - x₁), y₁ + t * (y₂ - y₁), z₁ + t * (z₂ - z₁))

        on_segment = isapprox(
            norm(point1 .- point) + norm(point .- point2), norm(point1 .- point2)
        )
        at_start = isapprox(norm(point1 .- point), 0)
        at_end = isapprox(norm(point .- point2), 0)
        between = on_segment && !at_start && !at_end

        if !only_between || between
            push!(out, point)
        end
    end

    return out
end

"""
    gs_gs_distance(gs1::InputGS, gs2::InputGS) -> Real
    gs_gs_distance(::Val{:ECEF}, gs1::InputGS, gs2::InputGS) -> Real
    gs_gs_distance(::Val{:lat_lon}, gs1::InputGS, gs2::InputGS) -> Real

Compute the distance [m] between two ground stations `gs1` and `gs2`.

The first method uses geodetic coordinates directly. The `:ECEF` variant first converts from
ECEF to geodetic. The `:lat_lon` variant is equivalent to the default method but with
explicit coordinate type specification.

The default method applies line-of-sight Euclidean distance if available, otherwise uses
the [haversine formula](https://en.wikipedia.org/wiki/Haversine_formula).

# Returns

  - `Real`: Distance [m].

# Throws

  - `ArgumentError`: If any coordinate is not finite.
"""
function gs_gs_distance(gs1::InputGS, gs2::InputGS)
    validate_coordinates(gs1)
    validate_coordinates(gs2)

    intersection_points = ellipsoid_line_intersection(
        SEMIMAJOR_RADIUS,
        SEMIMINOR_RADIUS,
        geodetic_to_ecef(gs1),
        geodetic_to_ecef(gs2);
        only_between=true,
    )
    both_on_surface =
        (length(gs1) == 2 || isapprox(gs1[3], 0; atol=5)) &&
            (length(gs2) == 2 || isapprox(gs2[3], 0; atol=5))

    if length(intersection_points) == 0 && !both_on_surface
        # use line of sight distance if there's line of sight between them
        return norm(geodetic_to_ecef(gs1) .- geodetic_to_ecef(gs2))
    end

    # otherwise, use haversine formula to compute distance along the surface of the Earth
    lat1, lon1 = gs1[1:2]
    lat2, lon2 = gs2[1:2]
    havθ = (1 - cos(lat2 - lat1) + cos(lat1) * cos(lat2) * (1 - cos(lon2 - lon1))) / 2
    θ = acos(1 - havθ * 2)
    return θ * (SEMIMAJOR_RADIUS + SEMIMINOR_RADIUS) / 2
end
function gs_gs_distance(::Val{:ECEF}, gs1::InputGS, gs2::InputGS)
    return gs_gs_distance(ecef_to_geodetic(gs1), ecef_to_geodetic(gs2))
end
function gs_gs_distance(::Val{:lat_lon}, gs1::InputGS, gs2::InputGS)
    return gs_gs_distance(gs1, gs2)
end

"""
    sat_sat_distance(
        sat_prop1::OrbitPropagatorSgp4, sat_prop2::OrbitPropagatorSgp4;
        time::Time=missing
    ) -> Union{Real,Nothing}
    sat_sat_distance(
        sat_tle1::TLE, sat_tle2::TLE;
        time::Time=missing
    ) -> Union{Real,Nothing}

Compute the distance [m] between two satellites, given by `sat_prop1` and `sat_prop2`
or `sat_tle1` and `sat_tle2`, at a specific `time` [Julian days or DateTime].

If `time` is not provided, the later epoch of the two satellites will be used.

See
[SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl)
for reference of `OrbitPropagatorSgp4` and
[SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)
for reference of `TLE`.

# Returns

  - `Union{Real,Nothing}`: Distance [m] or `nothing` if obstructed by Earth.
"""
function sat_sat_distance(
    sat_prop1::OrbitPropagatorSgp4, sat_prop2::OrbitPropagatorSgp4; time::Time=missing
)
    time = _get_time_JD(time, sat_prop1, sat_prop2)

    Δt1_s = _JD_to_seconds(time - sat_prop1.sgp4d.epoch)
    Δt2_s = _JD_to_seconds(time - sat_prop2.sgp4d.epoch)

    pos1 = Tuple(Propagators.propagate!(sat_prop1, Δt1_s)[1])
    pos2 = Tuple(Propagators.propagate!(sat_prop2, Δt2_s)[1])

    visibility_check = ellipsoid_line_intersection(
        SEMIMAJOR_RADIUS, SEMIMINOR_RADIUS, pos1, pos2
    )
    if length(visibility_check) > 0
        return nothing
    end

    return norm(pos1 .- pos2)
end
function sat_sat_distance(sat_tle1::TLE, sat_tle2::TLE; time::Time=missing)
    return sat_sat_distance(
        Propagators.init(Val(:SGP4), sat_tle1),
        Propagators.init(Val(:SGP4), sat_tle2);
        time=time,
    )
end

"""
    sun_position(time::Real) -> Point3D
    sun_position(time::DateTime) -> Point3D

Calculate the position of the Sun [m] at a specific `time` in Earth-centered Inertial (ECI)
coordinates.

# Arguments

  - `time` as one of:
      + `Real`: Julian days.
      + `DateTime`: DateTime object.

# Returns

  - `Point3D`: Position of the Sun [m] in ECI coordinates.
"""
function sun_position(time::Real)
    𝑇 = time - JD_J2000
    # https://astronomy.stackexchange.com/questions/28802/calculating-the-sun-s-position-in-eci
    # Calculate parameters
    # mean longitude, in radians
    𝐿 = (280.4606184 + ((36000.77005361 / 36525) * 𝑇)) * π / 180
    # mean anomaly, in radians
    𝑀 = (357.5277233 + ((35999.05034 / 36525) * 𝑇)) * π / 180
    # ecliptic longitude, in radians
    ℓ = 𝐿 + ((1.914666471 * sin(𝑀)) + (0.918994643 * sin(2 * 𝑀))) * π / 180
    # obliquity of ecliptic plane, in radians
    ε = (23.43929 - ((46.8093 / 3600) * (𝑇 / 36525))) * π / 180

    # Calculate unit directional vector in ECI coordinates
    û = [cos(ℓ), cos(ε) * sin(ℓ), sin(ε) * sin(ℓ)]

    # Calculate distance to sun and scale the unit vector
    # distance from Earth's center to Sun's center in astronomical units (AU)
    𝑅ₐ = 1.000140612 - (0.016708617 * cos(𝑀)) - (0.000139589 * cos(2 * 𝑀))
    # center-to-center distance from Earth to Sun in meters
    𝑅ₘ = 𝑅ₐ * 149597870700

    return 𝑅ₘ .* û # distance to sun in meters
end
function sun_position(time::DateTime)
    return sun_position(datetime2julian(time))
end

"""
    light_at_point(point_eci::Point3D, time::Union{Real, DateTime}) -> LightCondition

Determine the lighting condition at the Earth-Centered Inertial (ECI) coordinates given by
`point_eci` [m] at the specified `time` [Julian days or DateTime].

# Returns

  - `LightCondition`: Sunlight, penumbra, or umbra.

# Throws

  - `ArgumentError`: If any of the ECI coordinates are not finite.
"""
function light_at_point(point_eci::Point3D, time::Union{Real,DateTime})
    point_eci = _to_float64(point_eci)
    validate_coordinates(point_eci)

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