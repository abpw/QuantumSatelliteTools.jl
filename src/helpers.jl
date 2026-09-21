using Dates: DateTime, datetime2julian

using SatelliteToolboxBase: OrbitStateVector
using SatelliteToolboxPropagators: OrbitPropagatorSgp4

using ..constants: SECONDS_PER_DAY, LAND_SEA_MASK
using ..types: Point3D, InputGS, GS, Node, Time, satellite, ground_station

############################################################################################
#                                    Arithmetic Helpers                                    #
############################################################################################

"""
    norm(point::Point3D) -> Float64

Compute the Euclidean norm of a 3D point.
"""
function norm(point::Point3D)
    return √sum(point .^ 2)
end

"""
    dB_to_prob(loss_dB::Real) -> Float64

Convert a loss value [dB] to [prob].
"""
function dB_to_prob(loss_dB::Real)
    return 1 - 10^(loss_dB / -10)
end

"""
    prob_to_dB(loss_prob::Real) -> Float64

Convert loss value [prob] to [dB].
"""
function prob_to_dB(loss_prob::Real)
    return -10 * log10(1 - loss_prob)
end


############################################################################################
#                                 Type Conversion Helpers                                  #
############################################################################################

"""
    _to_float64(coords::InputGS) -> GS
    _to_float64(lat::Real, lon::Real) -> GS
    _to_float64(lat::Real, lon::Real, height::Real) -> GS

Convert `coords` or `(lat, lon [, height])` to Float64 precision.
"""
function _to_float64(coords::InputGS)
    return Float64.(coords)
end
function _to_float64(lat::Real, lon::Real)
    return (Float64(lat), Float64(lon))
end
function _to_float64(lat::Real, lon::Real, height::Real)
    return (Float64(lat), Float64(lon), Float64(height))
end

"""
    _parse_lat_lon(lat::AbstractString, lon::AbstractString) -> GS

Parse `lat` and `lon` [deg] to a ground station tuple [rad].
"""
function _parse_lat_lon(lat_str::AbstractString, lng_str::AbstractString)
    return (deg2rad(parse(Float64, lat_str)), deg2rad(parse(Float64, lng_str)))
end

"""
    _JD_to_seconds(jd::Real) -> Float64

Convert Julian days to seconds.
"""
function _JD_to_seconds(jd::Real)
    return jd * SECONDS_PER_DAY
end

"""
    _seconds_to_JD(seconds::Real) -> Float64

Convert seconds to Julian days.
"""
function _seconds_to_JD(seconds::Real)
    return seconds / SECONDS_PER_DAY
end

############################################################################################
#                                    Validation Helpers                                    #
############################################################################################

"""
    validate_non_negative_finite(name::AbstractString, value::Real) -> Nothing

Validate that `value` is non-negative and finite.

# Throws
  - `ArgumentError`: If `value` is negative or non-finite.
"""
function validate_non_negative_finite(name::AbstractString, value::Real)
    if value < 0 || !isfinite(value)
        throw(ArgumentError("$name must be finite and non-negative, got $value"))
    end
end

"""
    _validate_timeout(timeout::Union{Integer, Float64}) -> Nothing

Validate that `timeout` is non-negative and not NaN.

# Throws
  - `ArgumentError`: If `timeout` is negative or NaN.
"""
function _validate_timeout(timeout::Real)
    if timeout < 0 || isnan(timeout)
        throw(ArgumentError("timeout must be non-negative and not NaN, got $timeout"))
    end
end

"""
    validate_coordinates(lat::Real, lon::Real) -> Nothing
    validate_coordinates(lat::Real, lon::Real, height::Real) -> Nothing
    validate_coordinates(coords::InputGS) -> Nothing
    
Validate that the all values in `(lat, lon[, height])` or `coords` are finite numbers.

# Throws

  - `ArgumentError` if any coordinate is infinite or NaN.
"""
function validate_coordinates(lat::Real, lon::Real)
    if !isfinite(lat) || !isfinite(lon)
        throw(ArgumentError("Coordinates must be finite, got ($(lat), $(lon))."))
    end
end
function validate_coordinates(lat::Real, lon::Real, height::Real)
    if !isfinite(lat) || !isfinite(lon) || !isfinite(height)
        throw(ArgumentError("Coordinates must be finite, got ($(lat), $(lon), $(height))."))
    end
end
function validate_coordinates(coords::InputGS)
    validate_coordinates(coords...)
end

"""
    _to_JD(time::Time) -> Float64

Convert a `time` value to a Float64 Julian day.
"""
function _to_JD(time::Time)
    return typeof(time) == DateTime ? datetime2julian(time) : time
end

"""
    _get_time_JD(time::Time, sat_sv::OrbitStateVector) -> Float64
    _get_time_JD(time::Time, sat::OrbitPropagatorSgp4) -> Float64
    _get_time_JD(time::Time, sat1::OrbitPropagatorSgp4, sat2::OrbitPropagatorSgp4) -> Float64

Convert or populate `time` to a Float64 Julian day value.

# Arguments

  - `time::Time`: Time in Julian days or DateTime, or `missing` to use the satellite epoch.
  - `sat_sv::OrbitStateVector`: Satellite state vector (used if `time` is missing).
  - `sat::OrbitPropagatorSgp4`: Satellite propagator (used if `time` is missing).
  - `sat1::OrbitPropagatorSgp4`, `sat2::OrbitPropagatorSgp4`: Two satellite propagators
    (used if `time` is missing).

# Returns

  - `Float64`: Time in Julian days.
"""
function _get_time_JD(time::Time, sat_sv::OrbitStateVector)
    return time === missing ? sat_sv.t : _to_JD(time)
end
function _get_time_JD(time::Time, sat::OrbitPropagatorSgp4)
    return time === missing ? sat.sgp4d.epoch : _to_JD(time)
end
function _get_time_JD(time::Time, sat1::OrbitPropagatorSgp4, sat2::OrbitPropagatorSgp4)
    return time === missing ? max(sat1.sgp4d.epoch, sat2.sgp4d.epoch) : _to_JD(time)
end

############################################################################################
#                                    Condition Helpers                                     #
############################################################################################

"""
    is_satellite(node::Node) -> Bool

Check if a `node` is a satellite.
"""
function is_satellite(node::Node)
    return node.label.type == satellite
end

"""
    is_ground_station(node::Node) -> Bool

Check if a `node` is a ground station.
"""
function is_ground_station(node::Node)
    return node.label.type == ground_station
end

"""
    is_duplicate(new_gs::InputGS, gses::Vector; atol=1e-2) -> Bool

Check if `new_gs` matches an existing station in `gses` (within `atol` tolerance).

# Returns
  - `Bool`: `true` if `new_gs` matches an existing station, `false` otherwise.
"""
function is_duplicate(new_gs::InputGS, gses::Vector; atol=1e-2)
    new_coords = _to_float64(new_gs)
    return any(gs -> all(isapprox.(new_coords, _to_float64(gs); atol)), gses)
end

"""
    is_outside_min_distance(
        lat::Real, lon::Real, gses::Vector, min_distance::Real
    ) -> Bool

Check if `(lat, lon)` [rad] is at least `min_distance` [km] away from all `gses`.

# Returns
  - `Bool`: `true` if the location is sufficiently far from all stations.

# Throws
  - `ArgumentError`: If inputs are invalid (non-finite, out of range).
"""
function is_outside_min_distance(lat::Real, lon::Real, gses::Vector, min_distance::Real)
    _validate_non_negative_finite("min_distance", min_distance)

    coords = _to_float64(lat, lon)
    validate_coordinates(coords)

    gses = _to_float64.(gses)
    validate_coordinates.(gses)

    return all(gs -> gs_gs_distance((lat, lon), gs) >= min_distance * 1000, gses)
end

"""
    is_land(lat::Real, lon::Real) -> Bool
    is_land(coords::GS) -> Bool

Determine if geographic coordinates `(lat, lon)` or `coords` [rad] fall on land.

# Returns
  - `Bool`: `true` if coordinates are on land, `false` otherwise.

# Throws
  - `ArgumentError`: If coordinates are invalid (non-finite).
"""
function is_land(lat::Real, lon::Real)
    lat, lon = _to_float64(lat, lon)
    validate_coordinates((lat, lon))

    # convert radians to mask indices
    lat_idx = round(Int, (lat + π/2) / π * (length(LAND_SEA_MASK[2]) - 1) + 1)
    lon_idx = round(Int, (lon + π) / 2π * (length(LAND_SEA_MASK[1]) - 1) + 1)

    return LAND_SEA_MASK[3][lon_idx, lat_idx] == 1
end
function is_land(coords::GS)
    return is_land(coords[1], coords[2])
end