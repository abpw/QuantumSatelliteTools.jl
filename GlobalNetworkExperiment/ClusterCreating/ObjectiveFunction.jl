# PASSED TESTS

using NearestNeighbors
include(raw"C:\\Users\\alexb\\ProgrammingProjects\\URV2025\\QSATJulia\\QuantumSatelliteTools.jl\\GlobalNetworkExperiment\\ClusterCreating\\DistanceCalc.jl")

dist_km(city1, city2) = dist(city1, city2) / 1000.0
MAX_DIST_KM = 450.0

function objective_function(dist::Float64, Cities::Vector{NTuple{2,Float64}}) 
    return 1- dist/(1+graph_diameter_km(Cities).diameter_km)
end


# --- weighted adjacency: neighbors with edge length in km ---
function build_weighted_adjacency(coords::Vector{NTuple{2,Float64}})
    n = length(coords)
    adjMat = [Vector{Tuple{Int,Float64}}() for _ in 1:n]
    for i in 1:n-1, j in i+1:n
        dkm = dist_km(coords[i], coords[j])
        if dkm <= MAX_DIST_KM
            push!(adjMat[i], (j, dkm))
            push!(adjMat[j], (i, dkm))
        end
    end
    return adjMat
end


function dijkstra(adjw::Vector{Vector{Tuple{Int,Float64}}}, s::Int)
    n = length(adjw)
    dist = fill(Inf, n)
    seen = falses(n)
    dist[s] = 0.0
    for _ in 1:n
        u = 0; best = Inf
        @inbounds for v in 1:n
            if !seen[v] && dist[v] < best
                best = dist[v]; u = v
            end
        end
        (u == 0 || best == Inf) && break 
        seen[u] = true
        du = dist[u]
        @inbounds for (w, wlen) in adjw[u]
            newd = du + wlen
            if newd < dist[w]
                dist[w] = newd
            end
        end
    end
    return dist
end

function graph_diameter_km(coords::Vector{NTuple{2,Float64}})
    adjw = build_weighted_adjacency(coords)
    n = length(adjw)

    best_d_km = 0.0
    best_pair = (0,0)
    connected = true

    for s in 1:n
        dist = dijkstra(adjw, s)
        if any(isinf, dist)
            connected = false
        end
        for t in 1:n
            if isfinite(dist[t]) && dist[t] > best_d_km
                best_d_km = dist[t]
                best_pair = (s, t)
            end
        end
    end

    return (; diameter_km = best_d_km,
            pair = best_pair,
            adjacency_weighted = adjw,
            is_connected = connected)
end