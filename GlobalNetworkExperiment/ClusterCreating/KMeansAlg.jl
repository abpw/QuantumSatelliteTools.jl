using Distances
using StaticArrays

const DEFAULT_EARTH_RADIUS_M = 6_371_000.0
@inline deg2rad(x::Float64) = x * π/ 180

# Convert (lat, lon) in degrees to a 2 x n matrix in radians (rows = lat, lon; cols = points)
function _to_rad_matrix(pts::Vector{NTuple{2,Float64}})
    isempty(pts) && return Array{Float64}(undef, 2, 0)
    return hcat((SVector(deg2rad(p[1]), deg2rad(p[2])) for p in pts)...)
end

# Group items by 1...k assignments; returns Vector{Vector{T}}
function _bucket_by_assignments(items, assignments::Vector{Int}, k::Int)
    @assert length(items) == length(assignments)
    @assert all(1 .<= assignments .<= k) "Assignments must be in 1..$k"
    buckets = [eltype(items)[] for _ in 1:k]
    @inbounds for (x, c) in zip(items, assignments)
        push!(buckets[c], x)
    end
    return buckets
end

# Optional: find each head's index in Cities (exact tuple match).
# Returns Vector{Union{Int,Missing}} aligned with `cluster_heads`.
function _find_head_indices_in_cities(Cities::Vector{NTuple{2,Float64}}, cluster_heads::Vector{NTuple{2,Float64}})
    idxmap = Dict{NTuple{2,Float64},Int}()
    @inbounds for (i, c) in pairs(Cities)
        # keep the first index for duplicates; change to `idxmap[c] = i` to keep last
        get!(idxmap, c, i)
    end
    [get(idxmap, h, missing) for h in cluster_heads]
end



"""
    clustering_fixed_heads(
        Cities::Vector{NTuple{2,Float64}},
        cluster_heads::Vector{NTuple{2,Float64}};
        earth_radius_m::Float64 = DEFAULT_EARTH_RADIUS_M
    ) -> NamedTuple

Assign each city to the **nearest fixed head** via Haversine distance on a sphere of
radius `earth_radius_m`. Heads never move.

Returns a NamedTuple:
- `assignments::Vector{Int}`           : length `n`, labels in `1:k` (nearest head per city)
- `clusters::Vector{Vector{NTuple{2,Float64}}}` : `k` buckets of cities by assignment
- `head_coords::Vector{NTuple{2,Float64}}`     : the heads you passed in (fixed)
- `head_indices_in_cities::Vector{Union{Int,Missing}}` : index of each head in `Cities` if found
- `distances_km::Matrix{Float64}`      : `kxn` matrix of head→city distances in **km**
"""
function clustering_fixed_heads(Cities::Vector{NTuple{2,Float64}}, cluster_heads::Vector{NTuple{2,Float64}}; earth_radius_m::Float64=DEFAULT_EARTH_RADIUS_M)
    Crad = _to_rad_matrix(Cities)        
    Hrad = _to_rad_matrix(cluster_heads)  
    Dm = pairwise(Haversine(earth_radius_m), Hrad, Crad; dims=2)  # k x n
    Dkm = Dm ./ 1000.0
    k, n = size(Dkm)
    assignments = Vector{Int}(undef, n)

    # Assign each city to the nearest head
    @inbounds for j in 1:n
        _, idx = findmin(@view Dkm[:, j])
        assignments[j] = idx
    end

    clusters = _bucket_by_assignments(Cities, assignments, k)
    head_idxs = _find_head_indices_in_cities(Cities, cluster_heads)

    return (
        assignments = assignments,
        clusters = clusters,
        head_coords = cluster_heads,
        head_indices_in_cities = head_idxs,
        distances_km = Dkm,
    )
end

# --- Convenience wrapper that returns only the buckets ---
"""
    cluster_cities_fixed_heads(
        Cities::Vector{NTuple{2,Float64}},
        cluster_heads::Vector{NTuple{2,Float64}};
        earth_radius_m::Float64 = DEFAULT_EARTH_RADIUS_M
    ) -> Vector{Vector{NTuple{2,Float64}}}

Just the clustered city buckets (one vector per head), with heads fixed.
"""
function cluster_cities_fixed_heads(Cities::Vector{NTuple{2,Float64}}, cluster_heads::Vector{NTuple{2,Float64}}; earth_radius_m::Float64 = DEFAULT_EARTH_RADIUS_M
)
    result = clustering_fixed_heads(Cities, cluster_heads; earth_radius_m=earth_radius_m)
    return result.clusters
end