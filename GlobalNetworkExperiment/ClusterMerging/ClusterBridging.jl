using NearestNeighbors
include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\DistanceCalc.jl")

dist_km(city1, city2) = dist(city1, city2) / 1000.0
MAX_DIST_KM = 450.0
"""
interpolate_great_circle(p1::Tuple, p2::Tuple, t::Real) -> (lat, lon)

Great-circle interpolation between two (lat, lon) points (degrees).
t in [0,1].
"""
function interpolate_great_circle(p1::NTuple{2,Float64}, p2::NTuple{2,Float64}, t::Float64)
    lat1 = deg2rad(Float64(p1[1])); lon1 = deg2rad(Float64(p1[2]))
    lat2 = deg2rad(Float64(p2[1])); lon2 = deg2rad(Float64(p2[2]))

    # 3D unit vectors on the sphere
    u = (cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1))
    v = (cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2))

    uv_dotted = clamp(u[1]*v[1] + u[2]*v[2] + u[3]*v[3], -1.0, 1.0)
    d = acos(uv_dotted)

    if d ≈ 0.0
        return (Float64(p1[1]), Float64(p1[2]))
    end

    s = sin(d)
    a = sin((1 - t) * d) / s
    b = sin(t * d) / s

    x = a*u[1] + b*v[1]
    y = a*u[2] + b*v[2]
    z = a*u[3] + b*v[3]

    norm = sqrt(x^2 + y^2 + z^2)
    x /= norm; y /= norm; z /= norm

    lat = atan(z, sqrt(x^2 + y^2))   
    lon = atan(y, x)                  

    return (rad2deg(lat), rad2deg(lon))
end

"""
create_single_bridge_chain(p1, p2; max_distance=400.0) -> Vector{Tuple{Float64,Float64}}

Return intermediate bridge points along the great circle so that
each hop is <= `max_distance` km. Returns `[]` if none needed.
"""
function create_single_bridge_chain(p1::Tuple, p2::Tuple; max_distance::Real=MAX_DIST_KM)
    total = dist_km(p1,p2)
    if total <= max_distance
        return Tuple{Float64,Float64}[]
    end

    num_segments = ceil(Int, total / max_distance)   # >= 2 if total > max
    nodes = Tuple{Float64,Float64}[]
    for i in 1:num_segments-1
        t = i / num_segments
        push!(nodes, interpolate_great_circle(p1, p2, t))
    end
    return nodes
end



# Convert (lat, lon) in degrees -> 3D unit vector
function latlon_to_unit(p::NTuple{2,Float64})
    lat = deg2rad(Float64(p[1])) 
    lon = deg2rad(Float64(p[2])) 
    return (cos(lat)*cos(lon),  # x
            cos(lat)*sin(lon),  # y
            sin(lat))         # z
end

# Build a KDTree from a vector of (lat,lon) by mapping to 3 x N matrix of unit vectors
function build_kdtree_latlon(points::Vector{NTuple{2,Float64}})
    X = Matrix{Float64}(undef, 3, length(points))
    for (j, p) in enumerate(points)
        x, y, z = latlon_to_unit(p)
        @inbounds X[1,j] = x; X[2,j] = y; X[3,j] = z
    end
    return KDTree(X)
end

"""
find_closest_points_between_clusters(c1, c2) -> (p_from_c1, q_from_c2)

"""
function find_closest_points_between_clusters(c1::Vector{<:NTuple{2,<:Real}}, c2::Vector{<:NTuple{2,<:Real}})
    # Build the tree on the smaller set for speed
    if length(c2) < length(c1)
        tree_points, query_points, swapped = c2, c1, true
    else
        tree_points, query_points, swapped = c1, c2, false
    end

    tree = build_kdtree_latlon(tree_points)

    best_idx_tree = 0
    best_idx_query = 0
    best_dist = Inf

    x = Vector{Float64}(undef, 3) 
    for (i, p) in enumerate(query_points)
        x[1], x[2], x[3] = latlon_to_unit(p)
        idxs, dists = knn(tree, x, 1, true)
        d = dists[1]
        if d < best_dist
            best_dist = d
            best_idx_tree = idxs[1]
            best_idx_query = i
        end
    end

    a = tree_points[best_idx_tree]
    b = query_points[best_idx_query]
    return swapped ? (b, a) : (a, b)
end


"""
create_bridge_chains(clusters; distance_threshold=400.0) -> Vector{NTuple{2,Float64}}
For each pair of clusters, connect the closest points with a chain of bridge nodes
so that spacing less than or equal to 'distance_threshold' km. Returns all bridge nodes.
"""
function create_bridge_chains(clusters::Vector{Vector{NTuple{2,Float64}}}; distance_threshold::Float64=MAX_DIST_KM)
    bridge_nodes = Tuple{Float64,Float64}[]
    for i in 1:length(clusters)-1
        for j in i+1:length(clusters)
            p1, p2 = find_closest_points_between_clusters(clusters[i], clusters[j])
            append!(bridge_nodes,
                    create_single_bridge_chain(p1, p2; max_distance=distance_threshold))
        end
    end
    return bridge_nodes
end











# using NearestNeighbors
# using Clustering
# function naive_bridge_all_cluster(clustered_cities::KmeansResult)
#     bundled_cities = bundle_cities(clustered_cities)
#     anon_GSes = []
#     for bundle1 in bundled_cities
#         for bundle2 in bundle_cities
#             if bundle2 != bundle1
#                 city1,city2, _ = closest_cross_pair(bundle1, bundle2)
#                 anon_GS = midpoint_geodesic(city1[1],city1[2],city2[1],city2[2])
#                 push!(anon_GSes,anon_GS) 
#         end
#     end
#     return anon_GSes
# end



# function bundle_cities(res::Clustering.KmeansResult)
#     k = length(res.counts)
#     # preallocate exact sizes so we avoid push! growth
#     idxs = [Vector{Int}(undef, res.counts[c]) for c in 1:k]
#     fillpos = ones(Int, k)
#     @inbounds for (i, c) in pairs(res.assignments)  # i = point index, c = cluster id
#         p = fillpos[c]
#         idxs[c][p] = i
#         fillpos[c] = p + 1
#     end
#     return idxs  # Vector{Vector{Int}}; idxs[c] are the indices in cluster c
# end


# # Convert Vector{Vector} to d x N column matrix (no fancy allocations)
# function to_colmat(V::Vector{Vector{Float64}})
#     d = length(V[1])
#     T = promote_type(map(eltype, V)...)    # common numeric type
#     M = Array{T}(undef, d, length(V))
#     for j in 1:length(V)
#         M[:, j] = V[j]
#     end
#     return M
# end

# """
# closest_cross_pair(A, B)
# A, B :: Vector{<:AbstractVector{<:Real}} where each inner vector is a point (same dimension).
# Returns (a, b, dist): the closest pair a in A, b in B and their Euclidean distance.
# """
# function closest_cross_pair(A::Vector{Vector{Float64}}, B::Vector{Vector{Float64}})
#     if length(A) <= length(B)
#         tree = KDTree(to_colmat(A))                  # d x n
#         idxs, dists = knn(tree, to_colmat(B), 1)     # (1 x m) indices/distances
#         j = argmin(@view dists[1, :])
#         ia = idxs[1, j]
#         return (A[ia], B[j], dists[1, j])
#     else
#         tree = KDTree(to_colmat(B))                  # d x m
#         idxs, dists = knn(tree, to_colmat(A), 1)     # (1 x n)
#         i = argmin(@view dists[1, :])
#         ib = idxs[1, i]
#         return (A[i], B[ib], dists[1, i])
#     end
# end


# # returns (lat_mid, lon_mid) in degrees
# function midpoint_geodesic(lat1, lon1, lat2, lon2)
#     d1, l1 = deg2rad(lat1), deg2rad(lon1)
#     d2, l2 = deg2rad(lat2), deg2rad(lon2)

#     # unit vectors
#     x1, y1, z1 = cos(d1)*cos(l1), cos(d1)*sin(l1), sin(d1)
#     x2, y2, z2 = cos(d2)*cos(l2), cos(d2)*sin(l2), sin(d2)

#     x, y, z = x1 + x2, y1 + y2, z1 + z2
#     r = sqrt(x^2 + y^2 + z^2)
#     if r < 1e-12
#         error("Midpoint undefined for antipodal points (infinitely many).")
#     end
#     x /= r; y /= r; z /= r

#     # back to lat/lon
#     dm = atan(z, sqrt(x^2 + y^2))     # = atan2(z, √(x²+y²))
#     lm = atan(y, x)                   # = atan2(y, x)
#     latm = rad2deg(dm)
#     lonm = mod(rad2deg(lm) + 540, 360) - 180  # wrap to [-180,180)
#     return (latm, lonm)
# end

# end
