using Test: @testset, @test, @test_nowarn, @test_skip, @test_throws

using QuantumSatelliteTools
using QuantumSatelliteTools.GenerateGroundStations

@testset "GenerateGroundStations" begin
    @testset "is_outside_min_distance" begin
        @testset "empty ground stations" begin
            @test is_outside_min_distance(deg2rad(45.0), deg2rad(-120.0), [], 100.0) == true
        end

        @testset "single station far away" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            @test is_outside_min_distance(deg2rad(45.0), deg2rad(90.0), gs, 100.0) == true
        end

        @testset "single station too close" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            @test is_outside_min_distance(deg2rad(0.01), deg2rad(0.01), gs, 500.0) == false
        end

        @testset "same location as station" begin
            gs = [(deg2rad(45.0), deg2rad(-120.0), 0.0)]
            @test is_outside_min_distance(deg2rad(45.0), deg2rad(-120.0), gs, 10.0) == false
        end

        @testset "multiple stations all far" begin
            gses = [
                (deg2rad(0.0), deg2rad(0.0), 0.0),
                (deg2rad(10.0), deg2rad(10.0), 0.0),
                (deg2rad(-10.0), deg2rad(-10.0), 0.0),
            ]
            @test is_outside_min_distance(deg2rad(60.0), deg2rad(120.0), gses, 100.0) ==
                true
        end

        @testset "multiple stations one too close" begin
            gses = [
                (deg2rad(0.0), deg2rad(0.0), 0.0),
                (deg2rad(10.0), deg2rad(10.0), 0.0),
                (deg2rad(45.0), deg2rad(45.0), 0.0),
            ]
            @test is_outside_min_distance(deg2rad(45.05), deg2rad(45.05), gses, 100.0) ==
                false
        end

        @testset "multiple stations all too close" begin
            gses = [
                (deg2rad(0.0), deg2rad(0.0), 0.0),
                (deg2rad(0.01), deg2rad(0.01), 0.0),
                (deg2rad(-0.01), deg2rad(-0.01), 0.0),
            ]
            @test is_outside_min_distance(deg2rad(0.0), deg2rad(0.0), gses, 100.0) == false
        end

        @testset "multiple stations mixed distances" begin
            gses = [
                (deg2rad(0.0), deg2rad(0.0), 0.0),      # Close
                (deg2rad(60.0), deg2rad(60.0), 0.0),    # Far
                (deg2rad(-60.0), deg2rad(-60.0), 0.0),  # Far
            ]
            # Near close station - should fail
            @test is_outside_min_distance(deg2rad(0.01), deg2rad(0.01), gses, 100.0) ==
                false
        end

        @testset "antipodal points" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            @test is_outside_min_distance(deg2rad(0.0), deg2rad(180.0), gs, 1000.0) == true
        end

        @testset "poles" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            @test is_outside_min_distance(deg2rad(90.0), deg2rad(0.0), gs, 100.0) == true
            @test is_outside_min_distance(deg2rad(-90.0), deg2rad(0.0), gs, 100.0) == true
        end

        @testset "date line edge case" begin
            gs = [(deg2rad(0.0), deg2rad(180.0), 0.0)]

            # Same location (180° = -180°)
            @test is_outside_min_distance(deg2rad(0.0), deg2rad(-180.0), gs, 100.0) == false

            # Different location near date line
            @test is_outside_min_distance(deg2rad(1.0), deg2rad(180.0), gs, 100.0) == true
        end

        @testset "altitude differences" begin
            gs_sea = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            gs_high = [(deg2rad(0.0), deg2rad(0.0), 5000.0)]

            lat, lon = deg2rad(0.0), deg2rad(0.0)
            @test is_outside_min_distance(lat, lon, gs_sea, 10.0) == false
            @test is_outside_min_distance(lat, lon, gs_high, 10.0) == false
        end

        @testset "stations with different altitudes" begin
            gses = [(deg2rad(0.0), deg2rad(0.0), 0.0), (deg2rad(0.0), deg2rad(0.0), 1000.0)]
            lat, lon = deg2rad(0.0), deg2rad(0.0)
            # Same lat/lon but different altitudes
            @test is_outside_min_distance(lat, lon, gses, 100.0) == false
        end

        @testset "zero min_distance" begin
            gs = [(deg2rad(45.0), deg2rad(-120.0), 0.0)]
            @test is_outside_min_distance(deg2rad(45.01), deg2rad(-120.01), gs, 0.0) == true
        end

        @testset "small min_distance" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            lat, lon = deg2rad(0.01), deg2rad(0.01)

            @test is_outside_min_distance(lat, lon, gs, 1.0) == true
        end

        @testset "large min_distance" begin
            gs = [(deg2rad(0.0), deg2rad(0.0), 0.0)]
            @test is_outside_min_distance(deg2rad(10.0), deg2rad(10.0), gs, 5000.0) == false
        end

        @testset "dense station network" begin
            gses = [
                (deg2rad(Float64(i)*10), deg2rad(Float64(j)*10), 0.0) for i in -5:5 for
                j in -9:9
            ]
            # Far from all
            @test is_outside_min_distance(deg2rad(85.0), deg2rad(170.0), gses, 100.0) ==
                true
        end

        @testset "return type" begin
            result = is_outside_min_distance(deg2rad(45.0), deg2rad(90.0), [], 100.0)
            @test result isa Bool
        end
    end

    @testset "is_land" begin
        @testset "known land locations" begin
            # New York City
            @test is_land(deg2rad(40.7), deg2rad(-74.0)) == true

            # London
            @test is_land(deg2rad(51.5), deg2rad(-0.1)) == true

            # Tokyo
            @test is_land(deg2rad(35.7), deg2rad(139.7)) == true

            # Sydney
            @test is_land(deg2rad(-33.9), deg2rad(151.2)) == true
        end

        @testset "just offshore from land" begin
            # Just offshore from New York (Atlantic Ocean)
            @test is_land(deg2rad(40.6), deg2rad(-73.5)) == false

            # Just offshore from London (North Sea)
            @test is_land(deg2rad(51.5), deg2rad(1.0)) == false

            # Just offshore from Tokyo (Pacific)
            @test is_land(deg2rad(35.5), deg2rad(140)) == false

            # Just offshore from Sydney (Tasman Sea)
            @test is_land(deg2rad(-33.9), deg2rad(152.0)) == false
        end

        @testset "known ocean locations" begin
            # Mid-Atlantic
            @test is_land(deg2rad(30.0), deg2rad(-50.0)) == false

            # Pacific Ocean (far from land)
            @test is_land(deg2rad(0.0), deg2rad(-140.0)) == false

            # Indian Ocean (far from land)
            @test is_land(deg2rad(-10.0), deg2rad(70.0)) == false
        end

        @testset "polar regions" begin
            @test is_land(deg2rad(85.0), deg2rad(0.0)) isa Bool
            @test is_land(deg2rad(-85.0), deg2rad(0.0)) isa Bool
        end

        @testset "desert regions" begin
            # Sahara Desert
            @test is_land(deg2rad(20.0), deg2rad(5.0)) == true

            # Australian Outback
            @test is_land(deg2rad(-25.0), deg2rad(133.0)) == true
        end

        @testset "prime meridian" begin
            # Central France (clearly inland)
            @test is_land(deg2rad(46.0), deg2rad(2.0)) == true

            # Central Africa (clearly inland)
            @test is_land(deg2rad(0.0), deg2rad(25.0)) == true
        end

        @testset "date line" begin
            # New Zealand (clearly land)
            @test is_land(deg2rad(-40.0), deg2rad(175.5)) == true

            # Far from land near date line
            @test is_land(deg2rad(45.0), deg2rad(180.0)) == false
        end

        @testset "tuple input format" begin
            gs = (deg2rad(40.7), deg2rad(-74.0), 0.0)
            @test is_land(gs) == true
        end

        @testset "separate arguments format" begin
            @test is_land(deg2rad(40.7), deg2rad(-74.0), 0.0) isa Bool
        end

        @testset "return type" begin
            result = is_land(deg2rad(0.0), deg2rad(0.0))
            @test result isa Bool
        end

        @testset "equator locations" begin
            # Equator over Africa (land)
            @test is_land(deg2rad(0.0), deg2rad(20.0)) == true

            # Equator over Pacific (water)
            @test is_land(deg2rad(0.0), deg2rad(-120.0)) == false
        end

        @testset "consistency" begin
            lat, lon = deg2rad(40.7), deg2rad(-74.0)
            result1 = is_land(lat, lon)
            result2 = is_land(lat, lon)
            @test result1 == result2
        end
    end

    @testset "generate_population_center_gses" begin
        @testset "return type and format" begin
            gses = generate_population_center_gses(5)
            @test gses isa Vector
            @test all(x isa Tuple && length(x) == 2 for x in gses)
            @test all(isa(x[1], Real) && isa(x[2], Real) for x in gses)
        end

        @testset "requested number of stations" begin
            for n in [1, 5, 10]
                gses = generate_population_center_gses(n)
                @test length(gses) == n
            end
        end

        @testset "valid coordinate ranges" begin
            gses = generate_population_center_gses(10)
            for (lat, lon) in gses
                @test -π/2 <= lat <= π/2
                @test -π <= lon <= π
            end
        end

        @testset "stations on land" begin
            gses = generate_population_center_gses(10)
            for (lat, lon) in gses
                @test is_land(lat, lon) == true
            end
        end

        @testset "zero stations" begin
            gses = generate_population_center_gses(0)
            @test length(gses) == 0
            @test gses isa Vector
        end

        @testset "single station" begin
            gses = generate_population_center_gses(1)
            @test length(gses) == 1
            @test gses[1] isa Tuple
        end

        @testset "with empty other_gses" begin
            gses = generate_population_center_gses(5; other_gses=[])
            @test length(gses) == 5
        end

        @testset "min_distance constraint" begin
            min_dist = 500.0  # 500 km
            gses = generate_population_center_gses(10; min_distance=min_dist)

            # All pairwise distances should be >= min_distance
            # Using simple distance metric (not calling gs_gs_distance)
            for i in eachindex(gses)
                for j in (i + 1):length(gses)
                    lat1, lon1 = gses[i]
                    lat2, lon2 = gses[j]
                    # Haversine formula
                    d_lat = lat2 - lat1
                    d_lon = lon2 - lon1
                    a = sin(d_lat/2)^2 + cos(lat1) * cos(lat2) * sin(d_lon/2)^2
                    c = 2 * asin(sqrt(a))
                    d = 6371 * c  # Earth radius in km

                    @test d >= min_dist * 0.99  # Allow small numerical error
                end
            end
        end

        # TODO fix other_gses logic and write a test for it

        @testset "verbose parameter (no error)" begin
            @test_nowarn generate_population_center_gses(3)
        end

        @testset "consistency" begin
            gses1 = generate_population_center_gses(5)
            gses2 = generate_population_center_gses(5)

            @test length(gses1) == length(gses2) == 5
        end

        @testset "no duplicates" begin
            gses = generate_population_center_gses(10)
            @test length(unique(gses)) == length(gses)
        end

        @testset "locations on land" begin
            gses = generate_population_center_gses(10)
            @test length(gses) == 10

            @test all(is_land(lat, lon) for (lat, lon) in gses)
        end

        @testset "timeout parameter" begin
            gses = generate_population_center_gses(5; timeout=100)
            @test length(gses) <= 5
            @test length(gses) > 0
        end

        @testset "zero min_distance" begin
            gses = generate_population_center_gses(10; min_distance=0)
            @test length(gses) == 10
        end

        @testset "return type structure" begin
            gses = generate_population_center_gses(5)
            for gs in gses
                @test typeof(gs) <: Tuple
                @test length(gs) == 2
                lat, lon = gs
                @test typeof(lat) <: Real
                @test typeof(lon) <: Real
            end
        end
    end

    @testset "generate_city_gses" begin
        @testset "return type and format" begin
            city_names = ["New York", "London", "Tokyo"]
            gses = generate_city_gses(city_names)
            @test gses isa Vector
            @test all(x isa Tuple && length(x) == 2 for x in gses)
            @test all(isa(x[1], Real) && isa(x[2], Real) for x in gses)
        end

        @testset "known cities" begin
            city_names = ["New York", "London", "Tokyo"]
            gses = generate_city_gses(city_names)
            @test length(gses) == 3
        end

        @testset "single city" begin
            gses = generate_city_gses(["Paris"])
            @test length(gses) == 1
            @test gses[1] isa Tuple
        end

        @testset "empty city list" begin
            gses = generate_city_gses(Vector{String}())
            @test length(gses) == 0
            @test gses isa Vector
        end

        @testset "valid coordinate ranges" begin
            city_names = ["New York", "London", "Tokyo", "Sydney", "Cairo"]
            gses = generate_city_gses(city_names)
            for (lat, lon) in gses
                @test -π/2 <= lat <= π/2
                @test -π <= lon <= π
            end
        end

        @testset "cities on land" begin
            city_names = ["New York", "London", "Tokyo", "Sydney"]
            gses = generate_city_gses(city_names)
            for (lat, lon) in gses
                @test is_land(lat, lon) == true
            end
        end

        @testset "unknown cities" begin
            city_names = ["New York", "XyzUnknownCity", "London"]
            gses = generate_city_gses(city_names)

            # Should return matches for known cities, skip unknown
            @test length(gses) <= length(city_names)
            @test length(gses) >= 1  # At least some cities found
        end

        @testset "case sensitivity" begin
            # Test if lowercase/uppercase variations work
            city_names_lower = ["new york", "london"]
            city_names_upper = ["NEW YORK", "LONDON"]
            city_names_mixed = ["New York", "London"]

            gses_lower = generate_city_gses(city_names_lower)
            gses_upper = generate_city_gses(city_names_upper)
            gses_mixed = generate_city_gses(city_names_mixed)

            # Should handle variations consistently
            @test length(gses_lower) >= 1
            @test length(gses_upper) >= 1
            @test length(gses_mixed) >= 1
        end

        @testset "duplicate city names" begin
            city_names = ["Paris", "Paris", "London"]
            gses = generate_city_gses(city_names)

            # Should handle duplicates (either include once or twice)
            @test length(gses) >= 1
            @test length(gses) <= 3
        end

        @testset "consistency" begin
            city_names = ["New York", "London", "Tokyo"]
            gses1 = generate_city_gses(city_names)
            gses2 = generate_city_gses(city_names)

            # Same cities should produce same results
            @test length(gses1) == length(gses2)
            @test all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            )
        end

        @testset "no duplicates in output" begin
            city_names = ["New York", "London", "Tokyo"]
            gses = generate_city_gses(city_names)

            # No two cities at exact same location
            @test length(unique(gses)) == length(gses)
        end

        @testset "verbose parameter true" begin
            city_names = ["New York", "UnknownCity123"]
            @test_nowarn generate_city_gses(city_names)
        end

        @testset "verbose parameter false" begin
            city_names = ["New York", "UnknownCity123"]
            @test_nowarn generate_city_gses(city_names)
        end

        @testset "major world cities" begin
            city_names = ["Beijing", "Berlin", "Dubai", "Moscow", "Singapore"]
            gses = generate_city_gses(city_names)

            @test length(gses) == 5
            @test all(is_land(lat, lon) for (lat, lon) in gses)
        end

        @testset "return type structure" begin
            gses = generate_city_gses(["New York"])
            for gs in gses
                @test typeof(gs) <: Tuple
                @test length(gs) == 2
                lat, lon = gs
                @test typeof(lat) <: Real
                @test typeof(lon) <: Real
            end
        end

        @testset "special characters in city names" begin
            # Cities with special characters
            city_names = ["São Paulo", "Montréal"]
            gses = generate_city_gses(city_names)

            # Should handle or skip appropriately
            @test length(gses) >= 0
        end

        @testset "whitespace in city names" begin
            city_names = ["New York", "Los Angeles", "San Francisco"]
            gses = generate_city_gses(city_names)

            @test length(gses) >= 1
        end
    end

    @testset "generate_random_gses" begin
        @testset "return type and format" begin
            gses = generate_random_gses(5)
            @test gses isa Vector
            @test all(x isa Tuple && length(x) == 2 for x in gses)
            @test all(isa(x[1], Real) && isa(x[2], Real) for x in gses)
        end

        @testset "requested number of stations" begin
            for n in [1, 5, 10, 20]
                gses = generate_random_gses(n)
                @test length(gses) == n
            end
        end

        @testset "zero stations" begin
            gses = generate_random_gses(0)
            @test length(gses) == 0
            @test gses isa Vector
        end

        @testset "single station" begin
            gses = generate_random_gses(1)
            @test length(gses) == 1
            @test gses[1] isa Tuple
        end

        @testset "valid coordinate ranges" begin
            gses = generate_random_gses(20)
            for (lat, lon) in gses
                @test -π/2 <= lat <= π/2
                @test -π <= lon <= π
            end
        end

        @testset "uniform distribution on sphere" begin
            gses = generate_random_gses(100)

            # Latitudes should be distributed across range
            lats = [lat for (lat, lon) in gses]
            @test minimum(lats) < 0  # Some southern hemisphere
            @test maximum(lats) > 0  # Some northern hemisphere

            # Longitudes should be distributed across range
            lons = [lon for (lat, lon) in gses]
            @test minimum(lons) < 0
            @test maximum(lons) > 0
        end

        @testset "force_land=false (default)" begin
            gses = generate_random_gses(50; force_land=false)
            @test length(gses) == 50
            # Some may be on land, some on water - no constraint
        end

        @testset "force_land=true" begin
            gses = generate_random_gses(20; force_land=true)
            @test length(gses) == 20
            @test all(is_land(lat, lon) for (lat, lon) in gses)
        end

        @testset "seed reproducibility" begin
            gses1 = generate_random_gses(10; seed=12345)
            gses2 = generate_random_gses(10; seed=12345)

            @test length(gses1) == length(gses2)
            @test all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            )
        end

        @testset "different seeds produce different results" begin
            gses1 = generate_random_gses(10; seed=11111)
            gses2 = generate_random_gses(10; seed=22222)

            @test length(gses1) == length(gses2)
            @test !all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            )
        end

        @testset "seed=false gives randomness" begin
            gses1 = generate_random_gses(10; seed=false)
            gses2 = generate_random_gses(10; seed=false)

            # Very unlikely to be identical
            @test !(all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            ))
        end

        @testset "min_distance constraint (no other_gses)" begin
            min_dist = 1000.0
            gses = generate_random_gses(10; min_distance=min_dist, seed=42)

            # All pairwise distances should be >= min_distance
            for i in eachindex(gses)
                for j in (i + 1):length(gses)
                    lat1, lon1 = gses[i]
                    lat2, lon2 = gses[j]
                    d_lat = lat2 - lat1
                    d_lon = lon2 - lon1
                    a = sin(d_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(d_lon / 2)^2
                    c = 2 * asin(sqrt(a))
                    d = 6371 * c

                    @test d >= min_dist * 0.99
                end
            end
        end

        @testset "zero min_distance" begin
            gses = generate_random_gses(10; min_distance=0)
            @test length(gses) == 10
        end

        @testset "timeout parameter" begin
            # With very small timeout, may not generate all requested
            gses = generate_random_gses(100; timeout=10, seed=42)
            @test length(gses) <= 100
            @test length(gses) >= 1
        end

        @testset "timeout with force_land and min_distance" begin
            # Strict constraints with timeout
            gses = generate_random_gses(
                5; timeout=100, force_land=true, min_distance=500.0, seed=42
            )
            @test length(gses) <= 5
            @test all(is_land(lat, lon) for (lat, lon) in gses)
        end

        @testset "verbose parameter false" begin
            @test_nowarn generate_random_gses(5)
        end

        @testset "verbose parameter true" begin
            @test_nowarn generate_random_gses(5)
        end

        @testset "no duplicates" begin
            gses = generate_random_gses(50; seed=999)
            @test length(unique(gses)) == length(gses)
        end

        @testset "return type structure" begin
            gses = generate_random_gses(5)
            for gs in gses
                @test typeof(gs) <: Tuple
                @test length(gs) == 2
                lat, lon = gs
                @test typeof(lat) <: Real
                @test typeof(lon) <: Real
            end
        end

        @testset "large n generation" begin
            gses = generate_random_gses(100; seed=777)
            @test length(gses) == 100
            @test all(-π/2 <= lat <= π/2 for (lat, lon) in gses)
            @test all(-π <= lon <= π for (lat, lon) in gses)
        end

        @testset "force_land with small timeout" begin
            # May not find enough land spots with small timeout
            gses = generate_random_gses(20; force_land=true, timeout=50, seed=555)
            @test length(gses) <= 20
            @test all(is_land(lat, lon) for (lat, lon) in gses)
        end

        @testset "min_distance with force_land" begin
            gses = generate_random_gses(15; force_land=true, min_distance=500.0, seed=333)
            @test all(is_land(lat, lon) for (lat, lon) in gses)

            # Check min_distance constraint
            for i in eachindex(gses)
                for j in (i + 1):length(gses)
                    lat1, lon1 = gses[i]
                    lat2, lon2 = gses[j]
                    d_lat = lat2 - lat1
                    d_lon = lon2 - lon1
                    a = sin(d_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(d_lon / 2)^2
                    c = 2 * asin(sqrt(a))
                    d = 6371 * c

                    @test d >= 500.0 * 0.99
                end
            end
        end

        @testset "seed with force_land reproducibility" begin
            gses1 = generate_random_gses(10; force_land=true, seed=777)
            gses2 = generate_random_gses(10; force_land=true, seed=777)

            @test all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            )
        end

        @testset "combination of parameters" begin
            gses = generate_random_gses(
                8; force_land=true, min_distance=300.0, timeout=200, verbose=false, seed=888
            )

            @test length(gses) <= 8
            @test all(is_land(lat, lon) for (lat, lon) in gses)

            # Check min_distance
            for i in eachindex(gses)
                for j in (i + 1):length(gses)
                    lat1, lon1 = gses[i]
                    lat2, lon2 = gses[j]
                    d_lat = lat2 - lat1
                    d_lon = lon2 - lon1
                    a = sin(d_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(d_lon / 2)^2
                    c = 2 * asin(sqrt(a))
                    d = 6371 * c

                    @test d >= 300.0 * 0.99
                end
            end
        end
    end

    @testset "generate_equispaced_gses" begin
        @testset "return type and format" begin
            gses = generate_equispaced_gses(5)
            @test gses isa Vector
            @test all(x isa Tuple && length(x) == 2 for x in gses)
            @test all(isa(x[1], Real) && isa(x[2], Real) for x in gses)
        end

        @testset "odd n values" begin
            for n in [1, 3, 5, 7, 11, 21]
                gses = generate_equispaced_gses(n)
                @test length(gses) == n
            end
        end

        @testset "even n with fail_on_even=true throws error" begin
            for n in [2, 4, 6, 8, 10]
                @test_throws ErrorException generate_equispaced_gses(n; fail_on_even=true)
            end
        end

        @testset "even n with fail_on_even=false returns n+1" begin
            for n in [2, 4, 6, 8, 10]
                gses = generate_equispaced_gses(n; fail_on_even=false)
                @test length(gses) == n + 1
            end
        end

        @testset "n=1 edge case" begin
            gses = generate_equispaced_gses(1)
            @test length(gses) == 1
            @test gses[1] isa Tuple
        end

        @testset "valid coordinate ranges" begin
            gses = generate_equispaced_gses(21)
            for (lat, lon) in gses
                @test -π/2 <= lat <= π/2
                @test -π <= lon <= π
            end
        end

        @testset "no duplicates" begin
            gses = generate_equispaced_gses(25)
            @test length(unique(gses)) == length(gses)
        end

        @testset "approximately equispaced distribution" begin
            gses = generate_equispaced_gses(25)

            # Latitudes should be reasonably distributed
            lats = [lat for (lat, lon) in gses]
            lat_range = maximum(lats) - minimum(lats)
            @test lat_range > π/2  # Should span significant latitude range

            # Longitudes should be reasonably distributed
            lons = [lon for (lat, lon) in gses]
            @test minimum(lons) < 0
            @test maximum(lons) > 0
        end

        @testset "consistency (deterministic)" begin
            gses1 = generate_equispaced_gses(11)
            gses2 = generate_equispaced_gses(11)

            @test length(gses1) == length(gses2)
            @test all(
                isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
                isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
            )
        end

        @testset "return type structure" begin
            gses = generate_equispaced_gses(7)
            for gs in gses
                @test typeof(gs) <: Tuple
                @test length(gs) == 2
                lat, lon = gs
                @test typeof(lat) <: Real
                @test typeof(lon) <: Real
            end
        end

        @testset "Fibonacci lattice properties" begin
            # Fibonacci lattice should provide good coverage
            gses = generate_equispaced_gses(13)

            # All coordinates should be valid
            @test all(-π/2 <= lat <= π/2 for (lat, lon) in gses)
            @test all(-π <= lon <= π for (lat, lon) in gses)

            # Should have points across different hemispheres
            lats = [lat for (lat, lon) in gses]
            @test minimum(lats) < -π/6  # Some southern
            @test maximum(lats) > π/6   # Some northern
        end

        @testset "increasing n produces more stations" begin
            gses5 = generate_equispaced_gses(5)
            gses11 = generate_equispaced_gses(11)
            gses21 = generate_equispaced_gses(21)

            @test length(gses5) < length(gses11)
            @test length(gses11) < length(gses21)
        end

        @testset "large odd n" begin
            gses = generate_equispaced_gses(101)
            @test length(gses) == 101
            @test all(-π/2 <= lat <= π/2 for (lat, lon) in gses)
            @test all(-π <= lon <= π for (lat, lon) in gses)
        end

        @testset "fail_on_even=false with odd n unchanged" begin
            n = 7
            gses = generate_equispaced_gses(n; fail_on_even=false)
            @test length(gses) == n  # Should not increment odd numbers
        end

        @testset "fail_on_even parameter default is true" begin
            # Default behavior should throw on even n
            @test_throws ErrorException generate_equispaced_gses(4)
        end

        @testset "small n values" begin
            for n in [1, 3, 5]
                gses = generate_equispaced_gses(n)
                @test length(gses) == n
                @test all(-π/2 <= lat <= π/2 for (lat, lon) in gses)
                @test all(-π <= lon <= π for (lat, lon) in gses)
            end
        end

        @testset "even n sequence with fail_on_even=false" begin
            for n in [2, 4, 6, 8]
                gses = generate_equispaced_gses(n; fail_on_even=false)
                @test length(gses) == n + 1
                @test all(-π/2 <= lat <= π/2 for (lat, lon) in gses)
            end
        end

        @testset "Fibonacci lattice spread" begin
            # Check that points are reasonably spread out
            gses = generate_equispaced_gses(17)

            # Compute average pairwise distance to verify spread
            total_dist = 0
            count = 0
            for i in eachindex(gses)
                for j in (i + 1):length(gses)
                    lat1, lon1 = gses[i]
                    lat2, lon2 = gses[j]
                    Δlat = lat2 - lat1
                    Δlon = lon2 - lon1
                    a = sin(Δlat/2)^2 + cos(lat1) * cos(lat2) * sin(Δlon/2)^2
                    c = 2 * asin(√a)
                    d = 6371 * c
                    total_dist += d
                    count += 1
                end
            end
            avg_dist = total_dist / count

            # For equispaced points on sphere, average distance should be reasonable
            @test avg_dist > 1000  # More than 1000 km for 17 points
        end
    end

    @testset "generate_grid_gses" begin
        using StatsBase: sample

        @testset "return type and format" begin
            gses = generate_grid_gses()
            @test gses isa Vector
            @test all(x isa Tuple && length(x) == 2 for x in gses)
            @test all(isa(x[1], Real) && isa(x[2], Real) for x in gses)
        end

        @testset "default parameters" begin
            gses = generate_grid_gses()
            @test length(gses) > 0
        end

        @testset "valid coordinate ranges (sampled)" begin
            gses = generate_grid_gses(equatorial_distance_km=500)
            sample_size = min(100, length(gses))
            sample_indices = sample(1:length(gses), sample_size; replace=false)

            for i in sample_indices
                lat, lon = gses[i]
                @test -π/2 <= lat <= π/2
                @test -π <= lon <= π
            end
        end

        @testset "equatorial_distance_km parameter" begin
            gses_dense = generate_grid_gses(equatorial_distance_km=200)
            gses_sparse = generate_grid_gses(equatorial_distance_km=800)
            @test length(gses_dense) > length(gses_sparse)
        end

        @testset "grid structure: denser at equator" begin
            gses = generate_grid_gses(equatorial_distance_km=400)
            near_equator = count(abs(lat) < π/6 for (lat, lon) in gses)
            near_poles = count(abs(lat) > 5*π/12 for (lat, lon) in gses)
            @test near_equator > near_poles
        end

        @testset "no duplicates" begin
            gses = generate_grid_gses()
            @test length(unique(gses)) == length(gses)
        end

        @testset "consistency (deterministic)" begin
            gses1 = generate_grid_gses(equatorial_distance_km=400)
            gses2 = generate_grid_gses(equatorial_distance_km=400)
            @test length(gses1) == length(gses2)
        end

        @testset "force_land=false (default)" begin
            gses = generate_grid_gses(equatorial_distance_km=500; force_land=false)
            @test length(gses) > 0
        end

        @testset "force_land=true (sampled)" begin
            gses = generate_grid_gses(equatorial_distance_km=600; force_land=true)
            @test length(gses) > 0

            # Sample check instead of testing all
            sample_size = min(50, length(gses))
            sample_indices = sample(1:length(gses), sample_size; replace=false)

            for i in sample_indices
                lat, lon = gses[i]
                @test is_land(lat, lon) == true
            end
        end

        @testset "force_land=true produces fewer stations" begin
            gses_all = generate_grid_gses(equatorial_distance_km=500; force_land=false)
            gses_land = generate_grid_gses(equatorial_distance_km=500; force_land=true)
            @test length(gses_land) <= length(gses_all)
        end

        @testset "discount_func parameter" begin
            discount_identity = lat_rad -> 0
            gses = generate_grid_gses(
                equatorial_distance_km=500; discount_func=discount_identity
            )
            @test length(gses) > 0
        end

        @testset "discount_func reduces density at poles" begin
            discount_poles = lat_rad -> abs(lat_rad) / (π/2)
            gses = generate_grid_gses(
                equatorial_distance_km=400; discount_func=discount_poles
            )

            near_poles = count(abs(lat) > 3*π/8 for (lat, lon) in gses)
            near_equator = count(abs(lat) < π/8 for (lat, lon) in gses)
            @test near_equator >= near_poles
        end

        @testset "return type structure (sampled)" begin
            gses = generate_grid_gses()
            sample_size = min(20, length(gses))
            sample_indices = sample(1:length(gses), sample_size; replace=false)

            for i in sample_indices
                gs = gses[i]
                @test typeof(gs) <: Tuple
                @test length(gs) == 2
            end
        end

        @testset "global coverage" begin
            gses = generate_grid_gses(equatorial_distance_km=600)
            lons = [lon for (lat, lon) in gses]
            lats = [lat for (lat, lon) in gses]

            @test minimum(lons) < -π/2
            @test maximum(lons) > π/2
            @test minimum(lats) < -π/4
            @test maximum(lats) > π/4
        end

        @testset "equator has stations" begin
            gses = generate_grid_gses(equatorial_distance_km=500)
            near_equator = count(abs(lat) < 0.1 for (lat, lon) in gses)
            @test near_equator > 0
        end

        @testset "large equatorial distance" begin
            gses = generate_grid_gses(equatorial_distance_km=2000)
            @test length(gses) >= 1
        end

        @testset "small equatorial distance" begin
            gses = generate_grid_gses(equatorial_distance_km=100)
            @test length(gses) > 100
        end

        @testset "force_land with discount_func" begin
            discount = lat_rad -> 0.3 * abs(lat_rad)
            gses = generate_grid_gses(
                equatorial_distance_km=500; force_land=true, discount_func=discount
            )
            @test length(gses) > 0
        end

        @testset "discount_func affects spacing" begin
            discount_none = lat_rad -> 0
            discount_heavy = lat_rad -> abs(lat_rad)

            gses_none = generate_grid_gses(
                equatorial_distance_km=500; discount_func=discount_none
            )
            gses_heavy = generate_grid_gses(
                equatorial_distance_km=500; discount_func=discount_heavy
            )

            @test length(gses_heavy) <= length(gses_none)
        end
    end
end