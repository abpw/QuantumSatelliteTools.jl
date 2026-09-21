using Random: Xoshiro, rand, seed!
using Unicode: normalize

using CSV: Rows

using ..constants: EQUATORIAL_CIRCUMFERENCE_KM, SIN_60, LAND_SEA_MASK
using ..helpers: _to_float64, _parse_lat_lon
using ..helpers: validate_non_negative_finite, _validate_timeout
using ..helpers: is_duplicate, is_land
using ..types: GS
using ..AstronomyGeometry: validate_coordinates, gs_gs_distance

""" City and population data from:
    @misc{
        Youderain_2021,
        title={World Cities Database},
        url={https://simplemaps.com/data/world-cities},
        journal={simplemaps},
        author={Youderain, Chris},
        year={2021},
        month={Jun}
    }
"""

"""
    CITY_DATA_FILE

Path to the CSV file containing city and population data.

Each row should contain at least the following fields:
  - `city_ascii`: City name in ASCII characters.
  - `lat`: Latitude [deg].
  - `lng`: Longitude [deg].
"""
const CITY_DATA_FILE = joinpath(@__DIR__, "../databases/worldcities.csv")

"""
    GS_GENERATION_RNG

Random number generator for ground station generation.
"""
const GS_GENERATION_RNG = Xoshiro(rand(0:(2^31-1)))

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
    validate_non_negative_finite("min_distance", min_distance)

    coords = _to_float64(lat, lon)
    validate_coordinates(coords)

    gses = _to_float64.(gses)
    validate_coordinates.(gses)

    return all(gs -> gs_gs_distance((lat, lon), gs) >= min_distance * 1000, gses)
end

"""
    generate_population_center_gses(
        n::Integer;
        other_gses::Vector=[],
        timeout::Union{Integer, Float64}=Inf,
        min_distance::Real=0,
        verbose::Bool=false,
    ) -> Vector{GS}
    generate_population_center_gses(
        ::Val{:rows}, n::Integer;
        other_gses::Vector=[],
        timeout::Union{Integer, Float64}=Inf,
        min_distance::Real=0,
        verbose::Bool=false,
    ) -> Vector

Generate ground stations at the `n` most populous cities.

The first method returns processed ground stations `[rad]`. The second method
(`:rows` variant) returns raw CSV rows (internal use).

# Keyword Arguments
  - `other_gses::Vector=[]`: Existing ground stations.
  - `timeout::Union{Integer, Float64}=Inf`: Max rows to scan.
  - `min_distance::Real=0`: Minimum inter-station distance [km].
  - `verbose::Bool=false`: Print messages if `true`.

# Returns
  - `Vector{GS}`: Ground station coordinates [rad] (includes `other_gses`).
  - `Vector`: CSV rows (`:rows` variant, does not include `other_gses`).

# Throws
  - `ArgumentError`: If parameters are invalid (non-finite, out of range).
"""
function generate_population_center_gses(
    n::Integer;
    other_gses::Vector=[],
    timeout::Union{Integer,Float64}=Inf,
    min_distance::Real=0,
    verbose::Bool=false,
)
    rows = generate_population_center_gses(
        Val(:rows), n;
        other_gses,
        timeout,
        min_distance,
        verbose
    )
    new_stations = [_parse_lat_lon(row.lat, row.lng) for row in rows]

    return vcat(new_stations, other_gses)
end
function generate_population_center_gses(
    ::Val{:rows},
    n::Integer;
    other_gses::Vector=[],
    timeout::Union{Integer,Float64}=Inf,
    min_distance::Real=0,
    verbose::Bool=false,
)
    validate_non_negative_finite("n", n)
    _validate_timeout(timeout)
    validate_non_negative_finite("min_distance", min_distance)

    gses::Vector{GS} = _to_float64.(other_gses)
    validate_coordinates.(gses)

    selected_rows = []
    count = 0

    for (i, row) in enumerate(Rows(CITY_DATA_FILE))
        if count >= n || i > timeout
            if i > timeout && verbose
                println("Timeout: generated $count of $n population centers")
            end
            return selected_rows
        end

        lat, lon = _parse_lat_lon(row.lat, row.lng)

        if min_distance == 0 || is_outside_min_distance(lat, lon, gses, min_distance)
            push!(gses, (lat, lon))
            push!(selected_rows, row)
            count += 1
        end
    end

    return selected_rows
end

"""
    generate_city_gses(
        city_names::Vector{String};
        other_gses::Vector=[],
        verbose::Bool=false
    ) -> Vector{GS}

Generate ground stations at cities specified by `city_names`.

# Keyword Arguments
  - `other_gses::Vector=[]`: Existing ground stations.
  - `verbose::Bool=false`: Print unmatched cities if `true`.

# Returns
  - `Vector{GS}`: Ground station coordinates [rad] (includes `other_gses`).

# Throws
  - `ArgumentError`: If `other_gses` contains invalid coordinates.
"""
function generate_city_gses(
    city_names::Vector{String}; other_gses::Vector=[], verbose::Bool=false
)
    gses::Vector{GS} = _to_float64.(other_gses)
    validate_coordinates.(gses)

    remaining = normalize.(lowercase.(city_names); stripmark=true)

    for row in Rows(CITY_DATA_FILE)
        if isempty(remaining)
            break
        end

        row_city_normalized = normalize(lowercase(row.city_ascii))
        if row_city_normalized ∈ remaining
            city_gs = _parse_lat_lon(row.lat, row.lng)

            if !is_duplicate(city_gs, gses)
                push!(gses, city_gs)
            elseif verbose
                println("Duplicate: $(row.city_ascii) already in station list, skipping")
            end

            filter!(name -> name ≠ row_city_normalized, remaining)
        end
    end

    if verbose && !isempty(remaining)
        println("Not found in database: $(join(remaining, ", "))")
    end

    return gses
end

"""
    generate_random_gses(
        n::Integer;
        other_gses::Vector=[],
        timeout::Union{Integer, Float64}=Inf,
        min_distance::Real=0,
        force_land::Bool=false,
        verbose::Bool=false,
        seed::Union{Integer, AbstractString, Bool}=false,
    ) -> Vector{GS}

Generate `n` randomly placed ground stations on a sphere.

# Keyword Arguments
  - `other_gses::Vector=[]`: Existing ground stations.
  - `timeout::Union{Integer, Float64}=Inf`: Maximum random draws.
  - `min_distance::Real=0`: Minimum inter-station distance [km].
  - `force_land::Bool=false`: Restrict stations to land if `true`.
  - `verbose::Bool=false`: Print messages if `true`.
  - `seed::Union{Integer, AbstractString, Bool}=false`: RNG seed (use `false` for random).

# Returns
  - `Vector{GS}`: Ground station coordinates [rad] (includes `other_gses`).

# Throws
  - `ArgumentError`: If parameters are invalid (non-finite, out of range).
"""
function generate_random_gses(
    n::Integer;
    other_gses::Vector=[],
    timeout::Union{Integer,Float64}=Inf,
    min_distance::Real=0,
    force_land::Bool=false,
    verbose::Bool=false,
    seed::Union{Integer,AbstractString,Bool}=false,
)
    validate_non_negative_finite("n", n)
    _validate_timeout(timeout)
    validate_non_negative_finite("min_distance", min_distance)

    seed = seed === false ? rand(1:(2^31-1)) : seed
    if verbose
        println("Generating $n random stations using seed $seed")
    end
    seed!(GS_GENERATION_RNG, seed)

    gses = _to_float64.(other_gses)
    validate_coordinates.(gses)

    generated = 0
    attempts = 0

    while generated < n && attempts < timeout
        attempts += 1

        # Uniform random sampling on sphere
        lat = acos(2 * rand(GS_GENERATION_RNG) - 1) - π/2
        lon = π * (2 * rand(GS_GENERATION_RNG) - 1)

        outside_min_distance = is_outside_min_distance(lat, lon, gses, min_distance)
        distance_ok = min_distance == 0 || outside_min_distance
        land_ok = !force_land || is_land(lat, lon)

        if distance_ok && land_ok
            push!(gses, _to_float64((lat, lon)))
            generated += 1
        end
    end

    if attempts >= timeout && verbose
        println("Timeout reached: $generated out of $n random ground stations generated")
    end

    return gses
end

"""
    generate_equispaced_gses(
        n::Integer; fail_on_even::Bool=true
    ) -> Vector{Tuple{Float64, Float64}}

Generate `n` approximately equispaced ground stations over the globe using a Fibonacci
lattice, an algorithm that only works correctly for odd `n`.

If `fail_on_even` is `true`, an error is thrown for even `n`. Otherwise, `n` is set to `n+1`.

# Returns

  - `Vector{Tuple{Float64, Float64}}`: `(latitude, longitude)` pairs [rad] for generated
    ground stations.

# Throws

  - `ErrorException`: If `n` is even and `fail_on_even` is `true`.
"""
function generate_equispaced_gses(n::Integer; fail_on_even::Bool=true)
    validate_non_negative_finite("n", n)

    if iseven(n) && fail_on_even
        throw(error("Fibonacci lattice algorithm only works for odd n, got n=$n"))
    elseif iseven(n)
        n += 1
    end

    gses = []

    N = n ÷ 2
    ϕ = (1 + √5) / 2

    for i in (-N):N
        lon_sign = sign(i) == 0 ? 1 : sign(i)

        lat = asin(2 * i / (2N + 1))
        lon = lon_sign * ((lon_sign * i) % ϕ) * 2π / ϕ

        if lon < -π
            lon += 2π
        end
        if lon > π
            lon -= 2π
        end

        push!(gses, (Float64(lat), Float64(lon)))
    end

    return gses
end

"""
    generate_grid_gses(;
        equatorial_distance_km::Integer=400,
        discount_func::Function=lat_rad->0,
        other_gses::Vector=[],
        force_land::Bool=false,
    ) -> Vector{Tuple{Float64, Float64}}

Generate ground stations over a spherical grid, denser near equator.

# Keyword Arguments

  - `equatorial_distance_km::Integer=400`: Distance between stations at the equator [km].
  - `discount_func::Function=lat_rad->0`: Function adjusting distance at latitudes.
  - `other_gses::Vector=[]`: Existing ground stations to append to.
  - `force_land::Bool=false`: Only include points on land.

# Returns

  - `Vector{Tuple{Float64, Float64}}`: `(lat, lon)` pairs [rad] for generated
    ground stations.
"""
function generate_grid_gses(;
    equatorial_distance_km::Integer=400,
    discount_func::Function=lat_rad->0,
    other_gses::Vector=[],
    force_land::Bool=false,
)
    validate_non_negative_finite("equatorial_distance_km", equatorial_distance_km)

    gses = copy(other_gses)

    num_gses_layer = ceil(EQUATORIAL_CIRCUMFERENCE_KM / equatorial_distance_km)
    latitudinal_layer_Δrad =
        SIN_60 * equatorial_distance_km / EQUATORIAL_CIRCUMFERENCE_KM * 2π
    latitude_rad = 0.0
    gs_layers = [
        (latitude_rad, 2π / num_gses_layer * gs - π) for gs in 0:(num_gses_layer-1)
    ]
    append!(gses, [gs for gs in gs_layers if !force_land || is_land(gs)])

    while latitudinal_layer_Δrad + latitude_rad < π/2
        latitude_rad = latitude_rad + latitudinal_layer_Δrad
        longitudinal_distance = equatorial_distance_km + discount_func(latitude_rad)
        num_gses_layer = ceil(
            EQUATORIAL_CIRCUMFERENCE_KM * cos(latitude_rad) / longitudinal_distance
        )
        latitudinal_layer_Δrad =
            SIN_60 * longitudinal_distance / EQUATORIAL_CIRCUMFERENCE_KM * 2π
        gs_layers = [
            (sign*latitude_rad, 2π/num_gses_layer*gs-π) for gs in 0:(num_gses_layer-1) for
            sign in (-1, 1)
        ]
        append!(gses, [gs for gs in gs_layers if !force_land || is_land(gs)])
    end

    return gses
end