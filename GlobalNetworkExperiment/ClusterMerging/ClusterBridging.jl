using NearestNeighbors
using Clustering
function naive_bridge_all_cluster(clustered_cities::KmeansResult)
    bundled_cities = bundle_cities(clustered_cities)
    anon_GSes = []
    for bundle1 in bundled_cities
        for bundle2 in bundle_cities
            if bundle2 != bundle1
                city1,city2, _ = closest_cross_pair(bundle1, bundle2)
                anon_GS = midpoint_geodesic(city1[1],city1[2],city2[1],city2[2])
                push!(anon_GSes,anon_GS) 
        end
    end
    return anon_GSes
end



function bundle_cities(res::Clustering.KmeansResult)
    k = length(res.counts)
    # preallocate exact sizes so we avoid push! growth
    idxs = [Vector{Int}(undef, res.counts[c]) for c in 1:k]
    fillpos = ones(Int, k)
    @inbounds for (i, c) in pairs(res.assignments)  # i = point index, c = cluster id
        p = fillpos[c]
        idxs[c][p] = i
        fillpos[c] = p + 1
    end
    return idxs  # Vector{Vector{Int}}; idxs[c] are the indices in cluster c
end


# Convert Vector{Vector} to d x N column matrix (no fancy allocations)
function to_colmat(V::Vector{Vector{Float64}})
    d = length(V[1])
    T = promote_type(map(eltype, V)...)    # common numeric type
    M = Array{T}(undef, d, length(V))
    for j in 1:length(V)
        M[:, j] = V[j]
    end
    return M
end

"""
closest_cross_pair(A, B)
A, B :: Vector{<:AbstractVector{<:Real}} where each inner vector is a point (same dimension).
Returns (a, b, dist): the closest pair a in A, b in B and their Euclidean distance.
"""
function closest_cross_pair(A::Vector{Vector{Float64}}, B::Vector{Vector{Float64}})
    if length(A) <= length(B)
        tree = KDTree(to_colmat(A))                  # d x n
        idxs, dists = knn(tree, to_colmat(B), 1)     # (1 x m) indices/distances
        j = argmin(@view dists[1, :])
        ia = idxs[1, j]
        return (A[ia], B[j], dists[1, j])
    else
        tree = KDTree(to_colmat(B))                  # d x m
        idxs, dists = knn(tree, to_colmat(A), 1)     # (1 x n)
        i = argmin(@view dists[1, :])
        ib = idxs[1, i]
        return (A[i], B[ib], dists[1, i])
    end
end


# returns (lat_mid, lon_mid) in degrees
function midpoint_geodesic(lat1, lon1, lat2, lon2)
    d1, l1 = deg2rad(lat1), deg2rad(lon1)
    d2, l2 = deg2rad(lat2), deg2rad(lon2)

    # unit vectors
    x1, y1, z1 = cos(d1)*cos(l1), cos(d1)*sin(l1), sin(d1)
    x2, y2, z2 = cos(d2)*cos(l2), cos(d2)*sin(l2), sin(d2)

    x, y, z = x1 + x2, y1 + y2, z1 + z2
    r = sqrt(x^2 + y^2 + z^2)
    if r < 1e-12
        error("Midpoint undefined for antipodal points (infinitely many).")
    end
    x /= r; y /= r; z /= r

    # back to lat/lon
    dm = atan(z, sqrt(x^2 + y^2))     # = atan2(z, √(x²+y²))
    lm = atan(y, x)                   # = atan2(y, x)
    latm = rad2deg(dm)
    lonm = mod(rad2deg(lm) + 540, 360) - 180  # wrap to [-180,180)
    return (latm, lonm)
end

end
