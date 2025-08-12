using NearestNeighbors
include("DistanceCalc.jl")

MAX_DIST::Int64 = 400 # in km

function objective_function(dist::Float64, Cities::Vector{NTuple{2,Float64}}) 
    return 1- dist/(1+diam(Cities))
end



# --- Graph diameter (in hops) under the 400 km edge rule ---
function graph_diameter(coords::Vector{NTuple{2,Float64}})
    adj = build_adjacency(coords)
    n = length(adj)
    best_d = -1
    best_pair = (0, 0)
    for s in 1:n
        dist = bfs(adj, s)
        for t in 1:n
            if dist[t] > best_d
                best_d = dist[t]
                best_pair = (s, t)
            end
        end
    end
    return (; diameter_hops = best_d, pair = best_pair, adjacency = adj)
end

# --- BFS from a single source; distances in hops (-1 = unreachable) ---
function bfs(adj::Vector{Vector{Int}}, s::Int)
    n = length(adj)
    dist = fill(-1, n)
    dist[s] = 0
    q = Vector{Int}(undef, n); head = 1; tail = 1
    q[tail] = s
    while head <= tail
        v = q[head]; head += 1
        dv1 = dist[v] + 1
        for u in adj[v]
            if dist[u] == -1
                dist[u] = dv1
                tail += 1
                q[tail] = u
            end
        end
    end
    return dist
end

# --- build adjacency: edge if distance < cutoff_km ---
function build_adjacency(coords::Vector{NTuple{2,Float64}})
    n = length(coords)
    adj = [Int64[] for _ in 1:n]
    for i in 1:n-1, j in i+1:n
        if have_R_sine(coords[i][1],coords[i][2], coords[j][1], coords[j][2]) < MAX_DIST
            push!(adj[i], j) 
            push!(adj[j], i)
        end
    end
    return adj
end