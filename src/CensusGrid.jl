const census_data_file_str = joinpath(@__DIR__, "../databases/hybrid_qn/census_nodes.csv")
const census_earth_radius_km = 6371.0

"""Load regional sites in radians: US census nodes or China city records (iso2=CN)."""
function census_ground_candidates(region::AbstractString="us")
    if region in ("china", "global")
        return [(deg2rad(parse(Float64, row.lat)), deg2rad(parse(Float64, row.lng)))
                for row in CSV.Rows(city_data_file_str) if region == "global" || (!ismissing(row.iso2) && row.iso2 == "CN")]
    end
    region == "us" || throw(ArgumentError("grid region must be us, china, or global"))
    return [(deg2rad(parse(Float64, row.lat)), deg2rad(parse(Float64, row.lon)))
            for row in CSV.Rows(census_data_file_str)]
end

"""
    generate_census_grid_gses([candidates]; equatorial_distance_km=400)

Construct an approximately square geographic grid and snap intersections to
existing census nodes. Rows start at the southernmost candidate with spacing
`d/R`; longitude columns are centered on the bounding-box midpoint with spacing
`d/(R*cos(latitude))`. The equatorial spacing control `d` is in km; R=6371 km.
Only intersections inside the bounding box are considered. Reject snaps farther
than `d/2` and deduplicate stations, preserving traversal order (first input wins
distance ties). Inputs/outputs are `(latitude, longitude)` radians.

Defaults use its contiguous-US census nodes. Bounds must not cross
the antimeridian or include poles. Station count follows from spacing, not a quota.
"""
generate_census_grid_gses(; equatorial_distance_km::Real=400, region::AbstractString="us") =
    generate_census_grid_gses(census_ground_candidates(region); equatorial_distance_km=equatorial_distance_km,
                             global_grid=region == "global")

"""Equator-anchored global rings, uniform longitude around each ring, and one target per pole."""
function global_census_grid_targets(d::Real)
    isfinite(d) && 0 < d <= π*census_earth_radius_km || throw(ArgumentError("invalid grid spacing"))
    delta = d/census_earth_radius_km
    targets = Tuple{Float64,Float64}[(-π/2, 0.0)]
    for k in -floor(Int, π/(2delta)):floor(Int, π/(2delta))
        lat = k*delta
        abs(lat) >= π/2 && continue
        count = max(1, round(Int, 2π*cos(lat)/delta))
        append!(targets, [(lat, -π + 2π*j/count) for j in 0:count-1])
    end
    push!(targets, (π/2,0.0))
    return targets
end

function generate_census_grid_gses(candidates;
                                  equatorial_distance_km::Real=400, global_grid::Bool=false)
    d = Float64(equatorial_distance_km)
    isfinite(d) && 0 < d <= π*census_earth_radius_km ||
        throw(ArgumentError("census grid distance must be finite and in (0, pi*6371] km"))
    isempty(candidates) && return Tuple{Float64,Float64}[]
    all(gs -> length(gs) == 2 && all(isfinite, gs) &&
        (global_grid ? -π/2 <= gs[1] <= π/2 : -π/2 < gs[1] < π/2) && -π <= gs[2] <= π, candidates) ||
        throw(ArgumentError("candidates must be finite (latitude, longitude) radians, excluding poles"))
    south, north = extrema(first.(candidates))
    west, east = extrema(last.(candidates))
    global_grid || east - west < π || throw(ArgumentError("census grid requires a regional, non-antimeridian bounding box"))
    middle = (west + east)/2
    delta = d/census_earth_radius_km
    # Chord distance is monotone in great-circle distance and stable at zero.
    xyz = [(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)) for (lat, lon) in candidates]
    max_chord_squared = 4*sin(delta/4)^2  # great-circle distance d/2
    selected = Tuple{Float64,Float64}[]
    seen = Set{Tuple{Float64,Float64}}()
    targets = global_grid ? global_census_grid_targets(d) : Tuple{Float64,Float64}[]
    for lat in (global_grid ? Float64[] : south:delta:north)
        lon_step = delta/cos(lat)
        first_column = ceil(Int, (west-middle)/lon_step)
        last_column = floor(Int, (east-middle)/lon_step)
        for column in first_column:last_column
            lon = middle + column*lon_step
            push!(targets, (lat,lon))
        end
    end
    for (lat,lon) in targets
            target = (cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
            best_index, best_distance = 0, Inf
            for (i, point) in enumerate(xyz)
                # No point outside this latitude band can pass the d/2 cutoff.
                global_grid && abs(candidates[i][1]-lat) > delta/2 && continue
                distance = sum((point[j]-target[j])^2 for j in 1:3)
                if distance < best_distance
                    best_index, best_distance = i, distance
                end
            end
            if best_distance <= max_chord_squared
                gs = (Float64(candidates[best_index][1]), Float64(candidates[best_index][2]))
                if gs ∉ seen
                    push!(selected, gs)
                    push!(seen, gs)
                end
            end
    end
    return selected
end
