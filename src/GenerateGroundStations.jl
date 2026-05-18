using GeoDatasets
using CSV
using Random: Xoshiro, rand, seed!
using ..AstronomyGeometry: gs_gs_distance, GS, equatorial_circumference_km, sin_60

""" City and population data from:
@misc{Youderain_2021, title={World Cities Database}, url={https://simplemaps.com/data/world-cities}, journal={simplemaps}, author={Youderain, Chris}, year={2021}, month={Jun}}"""

const city_data_file_str = joinpath(@__DIR__, "../databases/worldcities.csv")

"""
Check if a given (lat, lon) is at least `min_distance` away from all ground stations in the list.

# Arguments
- `lat::Float64`: Latitude of the new candidate location (in radians).
- `lon::Float64`: Longitude of the new candidate location (in radians).
- `gses::Vector`: List of existing ground station tuples.
- `min_distance::Int`: Minimum allowed distance in kilometers.

# Returns
- `Bool`: `true` if the new location is sufficiently far from all existing stations.
"""
function is_within_min_distance(lat::Float64, lon::Float64, gses::Vector, min_distance::Int)
    new_gs = (lat, lon)
    for gs ∈ gses
        if gs_gs_distance(new_gs, gs) < min_distance * 1000
            return false
        end
    end
    return true
end

const 🌎 = GeoDatasets.landseamask(resolution='f', grid=1.25)

"""
Determine whether the specified geographic coordinates fall on land.

# Arguments
- `lat::Number`: Latitude in radians.
- `lon::Number`: Longitude in radians.

# Returns
- `Bool`: `true` if the point lies on land, `false` otherwise.
"""
function is_land(lat::Number, lon::Number, args...)
    lat_index = round(Int, (lat + π/2) / π * (length(🌎[2]) - 1) + 1)
    lon_index = round(Int, (lon + π) / 2π * (length(🌎[1]) - 1) + 1)
    if 🌎[3][lon_index, lat_index] == 1
        return true
    end
    return false
end
"""
Wrapper method for `is_land` that accepts a ground station tuple (geodetic latitude in radians, geodetic longitude in radians[, heignt in meters]).
"""
function is_land((lat, lon, args...)::GS)
    return is_land(lat, lon)
end

"""
Generate ground stations at the most populous cities based on a database.

# Arguments
- `n::Int`: Number of population centers to select.

# Keyword Arguments
- `other_gses::Vector=[]`: Existing ground stations to preserve.
- `timeout::Union{Integer, Float64}=Inf`: Max number of rows to scan from file.
- `min_distance::Number=0`: Minimum distance in meters required between any two stations.
- `verbose::Bool=false`: Print messages if true.

# generate_population_center_gses returns
- `Vector`: A vector of `(geodetic latitude in radians, geodetic longitude in radians)` tuples for selected ground stations.

# _generate_population_center_gses returns
- `Vector`: A vector of CSV rows with all available information for selected ground stations.
"""
function _generate_population_center_gses(n::Int; other_gses::Vector=[], timeout::Union{Integer, Float64}=Inf, min_distance::Number=0, verbose::Bool=false)
    gses = copy(other_gses)
    n_cities = 0
    for (i, row) ∈ enumerate(CSV.Rows(city_data_file_str))
        if n_cities >= n
            return gses
        end
        if i > timeout
            if verbose
                println("Timeout reached: $n_cities out of $n population center ground stations generated")
            end
            return gses
        end
        lat = π/180*parse(Float64, row.lat)
        lon = π/180*parse(Float64, row.lng)
        if min_distance == 0 || is_within_min_distance(lat, lon, gses, min_distance)
            push!(gses, row)
            n_cities += 1
        end
    end
end

"""
Generate ground stations at the most populous cities based on a database.

# Arguments
- `n::Int`: Number of population centers to select.

# Keyword Arguments
- `other_gses::Vector=[]`: Existing ground stations to preserve.
- `timeout::Union{Integer, Float64}=Inf`: Max number of rows to scan from file.
- `min_distance::Number=0`: Minimum distance in meters required between any two stations.
- `verbose::Bool=false`: Print messages if true.

# generate_population_center_gses returns
- `Vector`: A vector of `(geodetic latitude in radians, geodetic longitude in radians)` tuples for selected ground stations.

# _generate_population_center_gses returns
- `Vector`: A vector of CSV rows with all available information for selected ground stations.
"""
function generate_population_center_gses(n::Int; other_gses::Vector=[], timeout::Union{Integer, Float64}=Inf, min_distance::Number=0, verbose::Bool=false)
    return [(π/180*parse(Float64, gs.lat), π/180*parse(Float64, gs.lng)) for gs in _generate_population_center_gses(n; other_gses=other_gses, timeout=timeout, min_distance=min_distance, verbose=verbose)]
end

"""
Generate ground stations at specified city names.

# Arguments
- `city_names::Vector{String}`: List of city names to look for.

# Keyword Arguments
- `other_gses::Vector=[]`: Existing ground stations to preserve.
- `verbose::Bool=false`: Print missing cities if true.

# Returns
- `Vector`: A vector of `(geodetic latitude in radians, geodetic longitude in radians)` tuples for matched cities.
"""
function generate_city_gses(city_names::Vector{String}; other_gses::Vector=[], verbose::Bool=false)
    gses = copy(other_gses)
    remaining = copy(city_names)
    for row ∈ CSV.Rows(city_data_file_str)
        row_city = row.city_ascii
        if row_city ∈ remaining
            push!(gses, (π/180*parse(Float64, row.lat), π/180*parse(Float64, row.lng)))
            filter!(x -> x ≠ row_city, remaining)
        end
        if length(remaining) == 0
            return gses
        end
    end
    if verbose
        println("The following cities were not in the database, so no ground stations were generated for them: $(join(remaining, ", "))")
    end
    return gses
end

const gs_generation_rng = Xoshiro(rand(0:2^31-1))

"""
Generate randomly placed ground stations using uniform sampling on a sphere.

# Arguments
- `n::Int`: Number of random stations to generate.

# Keyword Arguments
- `other_gses::Vector=[]`: Existing ground stations to preserve.
- `timeout::Union{Integer, Float64}=Inf`: Max number of random draws before giving up.
- `min_distance::Number=0`: Minimum separation between new and existing stations.
- `force_land::Bool=false`: Restrict stations to land if true.
- `verbose::Bool=false`: Print generation details if true.
- `seed=false`: RNG seed (use false to randomly generate one).

# Returns
- `Vector`: Generated `(geodetic latitude in radians, geodetic longitude in radians)` tuples.
"""
function generate_random_gses(n::Int; other_gses::Vector=[], timeout::Union{Integer, Float64}=Inf, min_distance::Number=0, force_land::Bool=false, verbose::Bool=false, seed=false)
    if seed == false
        seed = rand(0:2^31-1)
    end
    if verbose
        println("Generating random ground stations using seed $seed")
    end
    seed!(gs_generation_rng, seed)
    gses = copy(other_gses)
    num_new_gses = 0
    while num_new_gses < n && timeout > 0
        timeout -= 1
        lat = acos(2rand(gs_generation_rng) - 1) - π/2
        lon = π * (2rand(gs_generation_rng) - 1)
        if (min_distance == 0 || is_within_min_distance(lat, lon, gses, min_distance)) && (!force_land || is_land(lat, lon))
            num_new_gses += 1
            push!(gses, (lat, lon))
        end
    end
    if timeout == 0 && verbose
        println("Timeout reached: $num_new_gses out of $n random ground stations generated")
    end
    return gses
end

"""
Generate `n` approximately equispaced ground stations over the globe using a Fibonacci lattice, an algorithm that only works correctly for odd n.

# Arguments
- `n::Int`: Number of ground stations to generate.
- `fail_on_even::Bool=true`: If true, throws an error if `n` is odd. Otherwise, sets `n=n+1`.

# Returns
- `Vector`: A list of `(geodetic latitude in radians, geodetic longitude in radians)` tuples.
"""
function generate_equispaced_gses(n::Int; fail_on_even::Bool=true)
    if n%2 == 0 && fail_on_even
        throw("Fibonacci lattice algorithm only works for odd n, given n=$n")
    elseif n%2 == 0
        n += 1
    end
    gses = []
    N = n÷2
    ϕ = (1 + √5)/2
    for i ∈ -N:N
        sign = Int(i>=0)*2-1
        lat = asin(2*i/(2*N+1))
        lon = sign*((sign*i)%ϕ)*2π/ϕ
        if lon < -π
            lon += 2π
        end
        if lon > π
            lon -= 2π
        end
        push!(gses, (lat, lon))
    end
    return gses
end

"""
Generate ground stations over a spherical grid, denser near equator.

# Keyword Arguments
- `equatorial_distance_km::Int=400`: Distance between stations at the equator.
- `discount_func::Function=lat_rad->0`: Function adjusting distance at latitudes.
- `other_gses::Vector=[]`: Existing ground stations to append to.
- `force_land::Bool=false`: Only include points on land.

# Returns
- `Vector`: A list of `(geodetic latitude in radians, geodetic longitude in radians)` tuples.
"""
function generate_grid_gses(;equatorial_distance_km::Int=400, discount_func::Function=lat_rad->0, other_gses::Vector=[], force_land::Bool=false)
    gses = copy(other_gses)
    num_gses_layer = ceil(equatorial_circumference_km/equatorial_distance_km)
    latitudinal_layer_delta_rad = sin_60*equatorial_distance_km/equatorial_circumference_km*2π
    latitude_rad = 0
    gs_layers = [(latitude_rad, 2π/num_gses_layer*gs-π) for gs ∈ 0:num_gses_layer-1]
    append!(gses, [gs for gs ∈ gs_layers if !force_land || is_land(gs)])
    while latitudinal_layer_delta_rad + latitude_rad < π/2
        latitude_rad = latitude_rad + latitudinal_layer_delta_rad
        longitudinal_distance = equatorial_distance_km + discount_func(latitude_rad)
        num_gses_layer = ceil(equatorial_circumference_km * cos(latitude_rad) / longitudinal_distance)
        latitudinal_layer_delta_rad = sin_60*longitudinal_distance / equatorial_circumference_km * 2π
        gs_layers = [(sign*latitude_rad, 2π/num_gses_layer*gs-π) for gs ∈ 0:num_gses_layer-1 for sign ∈ (-1, 1)]
        append!(gses, [gs for gs ∈ gs_layers if !force_land || is_land(gs)])
    end
    return gses
end