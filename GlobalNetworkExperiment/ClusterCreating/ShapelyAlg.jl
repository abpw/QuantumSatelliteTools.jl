using LinearAlgebra
include("DistanceCalc.jl")
include("ObjectiveFunction.jl")
NUM_OF_CITIES::Int64=20

"""
Generates a set of cluster heads, which will generate the number of clusters we need, allowing us to use K-Means to cluster

input: city tuple

output: vector of cluster heads
"""
function shapelyAlg(Cities::Vector{NTuple{2,Float64}}, alpha::Float64)
    cluster_heads::Vector{NTuple{2,Float64}} = []
    shapely_values::Vector{Float64} = []
    for curr_city in Cities
        curr_city_shapely_val::Float64 = 0.0;
        for index_city in Cities
            if (index_city != curr_city)
                # Need to figure out how to pass in command line arg (see Driver)
                curr_city_shapely_val += objective_function(dist(curr_city, index_city),Cities)
            end
        end 
        curr_city_shapely_val *= 1/2
        push!(shapely_values, curr_city_shapely_val)
    end 

    cities_cpy = Cities
    while length(cities_cpy) > 0
        cluster_head = maximum(shapely_values)
        push!(cluster_heads, cluster_head)
        local_cities = [curr_city for curr_city in Cities if dist(cluster_head, curr_city)>= alpha]
        cities_cpy = setdiff(cities_cpy,local_cities)
    end

    return cluster_heads
end

