using Test: @testset

using QuantumSatelliteTools

@testset "QuantumSatelliteTools.jl" begin
    include("test_astronomy_geometry.jl")
    include("test_freespace_channels.jl")
    include("test_generate_ground_stations.jl")
end