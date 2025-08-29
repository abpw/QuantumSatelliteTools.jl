using Test

# Include the individual modules directly
include("../GlobalNetworkExperiment/ClusterCreating/DistanceCalc.jl")
include("../GlobalNetworkExperiment/ClusterCreating/ObjectiveFunction.jl")
include("../GlobalNetworkExperiment/ClusterCreating/ShapelyAlg.jl")
include("../GlobalNetworkExperiment/ClusterCreating/KMeansAlg.jl")
include("../GlobalNetworkExperiment/ClusterMerging/ClusterBridging.jl")

@testset "GlobalNetworkExperiment Tests" begin

    @testset "DistanceCalc Tests" begin
        # Test basic distance calculation
        city1 = (40.7128, -74.0060)  # New York
        city2 = (34.0522, -118.2437) # Los Angeles

        distance = dist(city1, city2)
        @test distance > 0
        @test distance < 5000000  # Should be less than 5000 km in meters

        # Test same city distance
        @test dist(city1, city1) ≈ 0.0 atol = 1e-10

        # Test haversine distance function directly
        lat1, lon1 = 0.0, 0.0
        lat2, lon2 = 0.0, 1.0
        d = have_R_sine_distance(lat1, lon1, lat2, lon2)
        @test d > 0

        # Test antipodal points (should be approximately π * R)
        lat1, lon1 = 0.0, 0.0
        lat2, lon2 = 0.0, 180.0
        d = have_R_sine_distance(lat1, lon1, lat2, lon2)
        @test abs(d - π * 6_371_000.0) < 1000  # Within 1km tolerance
    end

    @testset "ShapelyAlg Tests" begin
        # Test with empty cities
        empty_cities = NTuple{2,Float64}[]
        result = shapelyAlg(empty_cities)
        @test result == NTuple{2,Float64}[]

        # Test with single city
        single_city = [(40.7128, -74.0060)]
        result = shapelyAlg(single_city)
        @test length(result) == 1
        @test result[1] == single_city[1]

        # Test with two cities
        two_cities = [(40.7128, -74.0060), (34.0522, -118.2437)]
        result = shapelyAlg(two_cities)
        @test length(result) >= 1
        @test length(result) <= 2

        # Test with multiple cities
        cities = [
            (40.7128, -74.0060),  # New York
            (34.0522, -118.2437), # Los Angeles
            (41.8781, -87.6298),  # Chicago
            (29.7604, -95.3698),  # Houston
            (39.9526, -75.1652)   # Philadelphia
        ]
        result = shapelyAlg(cities)
        @test length(result) >= 1
        @test length(result) <= length(cities)

        # Test that all returned heads are in original cities
        for head in result
            @test head in cities
        end
    end

    @testset "ObjectiveFunction Tests" begin
        # Test objective function with various distances
        cities = [(40.7128, -74.0060), (34.0522, -118.2437)]

        # Test with zero distance
        result = objective_function(0.0, cities)
        @test isfinite(result)  # Should be finite

        # Test with small distance
        result = objective_function(10.0, cities)
        @test isfinite(result)  # Should be finite

        # Test with large distance
        result = objective_function(1000.0, cities)
        @test isfinite(result)  # Should be finite

        # Test that function returns reasonable values
        # The function can be negative for large distances relative to graph diameter
        @test isfinite(result)

        # Test graph diameter calculation
        diameter_info = graph_diameter_km(cities)
        @test haskey(diameter_info, :diameter_km)
        @test haskey(diameter_info, :pair)
        @test haskey(diameter_info, :adjacency_weighted)
        @test haskey(diameter_info, :is_connected)
        @test diameter_info.diameter_km >= 0

        # Test build_weighted_adjacency
        adj = build_weighted_adjacency(cities)
        @test length(adj) == length(cities)
        @test all(x -> x isa Vector{Tuple{Int,Float64}}, adj)

        # Test dijkstra algorithm
        if !isempty(adj)
            distances = dijkstra(adj, 1)
            @test length(distances) == length(cities)
            @test distances[1] == 0.0  # Distance to self should be 0
        end
    end

    @testset "KMeansAlg Tests" begin
        # Test clustering with fixed heads
        cities = [
            (40.7128, -74.0060),  # New York
            (34.0522, -118.2437), # Los Angeles
            (41.8781, -87.6298),  # Chicago
            (29.7604, -95.3698),  # Houston
            (39.9526, -75.1652)   # Philadelphia
        ]

        cluster_heads = [(40.7128, -74.0060), (34.0522, -118.2437)]

        # Test clustering_fixed_heads
        result = clustering_fixed_heads(cities, cluster_heads)

        @test length(result.assignments) == length(cities)
        @test all(1 .<= result.assignments .<= length(cluster_heads))
        @test length(result.clusters) == length(cluster_heads)
        @test result.head_coords == cluster_heads
        @test size(result.distances_km) == (length(cluster_heads), length(cities))

        # Test that all cities are assigned to exactly one cluster
        all_assigned_cities = vcat(result.clusters...)
        @test length(all_assigned_cities) == length(cities)
        @test all(city in all_assigned_cities for city in cities)

        # Test cluster_cities_fixed_heads convenience function
        clusters = cluster_cities_fixed_heads(cities, cluster_heads)
        @test length(clusters) == length(cluster_heads)
        @test sum(length.(clusters)) == length(cities)

        # Test with empty cities
        empty_result = clustering_fixed_heads(NTuple{2,Float64}[], cluster_heads)
        @test empty_result.assignments == Int[]
        @test length(empty_result.clusters) == length(cluster_heads)
        @test all(isempty, empty_result.clusters)
        @test empty_result.distances_km == Matrix{Float64}(undef, length(cluster_heads), 0)

        # Test with single city
        single_city = [(40.7128, -74.0060)]
        single_result = clustering_fixed_heads(single_city, cluster_heads)
        @test length(single_result.assignments) == 1
        @test length(single_result.clusters) == length(cluster_heads)
    end

    @testset "ClusterBridging Tests" begin
        # Test cluster bridging functionality
        cluster1 = [(40.7128, -74.0060), (39.9526, -75.1652)]  # NY, Philly
        cluster2 = [(34.0522, -118.2437), (29.7604, -95.3698)] # LA, Houston

        # Test find_closest_points_between_clusters
        p1, p2 = find_closest_points_between_clusters(cluster1, cluster2)
        @test p1 in cluster1
        @test p2 in cluster2

        # Test create_single_bridge_chain with short distance
        short_chain = create_single_bridge_chain((40.0, -74.0), (40.1, -74.1))
        @test isempty(short_chain)  # Should be empty for short distances

        # Test create_single_bridge_chain with long distance
        long_chain = create_single_bridge_chain((40.0, -74.0), (34.0, -118.0))
        @test !isempty(long_chain)  # Should have bridge points for long distances

        # Test interpolate_great_circle
        p1 = (40.0, -74.0)
        p2 = (34.0, -118.0)
        mid_point = interpolate_great_circle(p1, p2, 0.5)
        @test typeof(mid_point) == NTuple{2,Float64}
        @test 34.0 <= mid_point[1] <= 40.0  # Latitude should be between
        @test -118.0 <= mid_point[2] <= -74.0  # Longitude should be between

        # Test create_bridge_chains
        clusters = [cluster1, cluster2]
        bridge_nodes = create_bridge_chains(clusters)
        @test typeof(bridge_nodes) == Vector{NTuple{2,Float64}}

        # Test that we can calculate distances between clusters
        for city1 in cluster1
            for city2 in cluster2
                distance = dist(city1, city2)
                @test distance > 0
                @test distance < 5000000  # Reasonable upper bound
            end
        end
    end










end

println("All GlobalNetworkExperiment tests completed!")
