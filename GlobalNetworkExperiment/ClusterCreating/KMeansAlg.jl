using LinearAlgebra
using Clustering

"""
Uses Julia's k-means algorithm
"""
function city_kmeans(Cities::Vector{NTuple{2,Float64}}, K::Int64)
    processed_cities = process_cities(Cities)
    clustered_cities = kmeans(processed_cities, K)
    return clustered_cities
end

"""
Processes cities, which are in the form of an vector of tuples into a 2 x NUM_OF_CITIES matrix
where each col is a lat (first row)-lon(second row) pair
"""
function process_cities(Cities::Vector{NTuple::{2,Float64}})
    return [col[i] for i in 1:2, col in Cities] 
end



