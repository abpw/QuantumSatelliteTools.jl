using Test: @testset, @test, @test_nowarn, @test_skip, @test_throws
using Statistics: mean

using StatsBase: sample
using Suppressor: @capture_out

using QuantumSatelliteTools
using QuantumSatelliteTools.constants
using QuantumSatelliteTools.types
using QuantumSatelliteTools.helpers
using QuantumSatelliteTools.GenerateGroundStations

# Helpers

function haversine_for_tests(lat1, lon1, lat2, lon2)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
    return 6371 * 2 * asin(sqrt(a))
end

function check_min_distance_constraint(gses, min_dist)
    for i in eachindex(gses)
        for j in (i+1):length(gses)
            d = haversine_for_tests(gses[i]..., gses[j]...)
            if d < min_dist * 0.99
                return false
            end
        end
    end
    return true
end

function is_approx_in(point, collection; atol=1e-3)
    return any(
        isapprox(point[1], p[1]; atol) && isapprox(point[2], p[2]; atol) for p in collection
    )
end

function in_coord_range(gs::GS)
    return -π/2 <= gs[1] <= π/2 && -π <= gs[2] <= π
end

function same_gses(gses1, gses2)
    return all(
        isapprox(gses1[i][1], gses2[i][1]; atol=1e-10) &&
            isapprox(gses1[i][2], gses2[i][2]; atol=1e-10) for i in eachindex(gses1)
    )
end

const NYC = (deg2rad(40.7), deg2rad(-73.9))
const LONDON = (deg2rad(51.5), deg2rad(-0.1))
const TOKYO = (deg2rad(35.7), deg2rad(139.7))
const SYDNEY = (deg2rad(-33.9), deg2rad(151.2))
const ST_PAUL = (deg2rad(44.95), deg2rad(-93.09))

const InputError = Union{ArgumentError,DomainError,MethodError,TypeError}

@testset "GenerateGroundStations" begin
    @testset "is_outside_min_distance" begin
        @testset "basic functionality" begin
            @test is_outside_min_distance(π/4, -2π/3, [], 100) == true
            @test is_outside_min_distance(π/4, π/2, [(0, 0, 0)], 100) == true
            @test is_outside_min_distance(0, 0, [(1.7e-4, 1.7e-4, 0)], 500) == false
            @test is_outside_min_distance(π/4, -2π/3, [(π/4, -2π/3, 0)], 10) == false
        end

        @testset "multiple stations" begin
            gses = [(0, 0, 0), (π/18, π/18, 0.0), (-π/18, π/18, 0)]
            @test is_outside_min_distance(π/3, 2π/3, gses, 100) == true
            @test is_outside_min_distance(1.7e-4, 1.7e-4, gses, 100) == false
        end

        @testset "geographic edge cases" begin
            @test is_outside_min_distance(0, π, [(0, 0, 0)], 1000) == true    # Antipodal
            @test is_outside_min_distance(π/2, 0, [(0, 0, 0)], 100) == true   # North Pole
            @test is_outside_min_distance(0, -π, [(0, π, 0)], 100) == false   # Date line
        end

        @testset "distance thresholds" begin
            @test is_outside_min_distance(
                deg2rad(45.01), deg2rad(-120.01), [(π/4, -2π/3, 0)], 0
            ) == true
            @test is_outside_min_distance(1.7e-4, 1.7e-4, [(0, 0, 0)], 1) == true
            @test is_outside_min_distance(π/18, π/18, [(0, 0, 0)], 5000) == false
        end

        @testset "altitude variations" begin
            @test is_outside_min_distance(0, 0, [(0, 0, 0)], 10) == false
            @test is_outside_min_distance(0, 0, [(0, 0, 5000)], 10) == false
            @test is_outside_min_distance(0, 0, [(0, 0, 0), (0, 0, 1000)], 100) == false
        end

        @testset "dense station network" begin
            gses = [(deg2rad(i*10), deg2rad(i*10), 0) for i in -5:5 for j in -9:9]
            @test is_outside_min_distance(deg2rad(85), deg2rad(170), gses, 100) == true
        end

        @testset "return type" begin
            @test is_outside_min_distance(π/4, π/2, [], 100) isa Bool
        end

        @testset "consistency" begin
            lat, lon = deg2rad(40.7), deg2rad(-74.0)
            gses = [(deg2rad(51.5), deg2rad(-0.1), 0), (deg2rad(35.7), deg2rad(139.7), 0)]
            result1 = is_outside_min_distance(lat, lon, gses, 100)
            result2 = is_outside_min_distance(lat, lon, gses, 100)
            @test result1 == result2
        end

        @testset "error handling" begin
            for val in [-Inf, NaN, Inf]
                @test_throws InputError is_outside_min_distance(val, 0, [], 100)
                @test_throws InputError is_outside_min_distance(0, val, [], 100)
                @test_throws InputError is_outside_min_distance(0, 0, [], val)
                for gs in [(val, 0, 0), (0, val, 0), (0, 0, val)]
                    @test_throws InputError is_outside_min_distance(0, 0, [gs], 100)
                end
            end
            @test_throws InputError is_outside_min_distance(0, 0, [], -100)
        end
    end

    @testset "is_land" begin
        @testset "known locations" begin
            # Land
            @test is_land(NYC) == true
            @test is_land(LONDON) == true
            @test is_land(TOKYO) == true
            @test is_land(SYDNEY) == true
            # Ocean (just offshore)
            @test is_land(deg2rad(40.6), deg2rad(-73.5)) == false  # Near NYC
            @test is_land(deg2rad(51.5), deg2rad(1)) == false      # Near London
            @test is_land(deg2rad(35.5), deg2rad(140)) == false    # Near Tokyo
            @test is_land(deg2rad(-33.9), deg2rad(152)) == false   # Near Sydney

            # Ocean
            @test is_land(π/6, deg2rad(-50)) == false              # Mid-Atlantic
            @test is_land(0, deg2rad(-140)) == false               # Pacific Ocean
            @test is_land(π/18, deg2rad(70)) == false              # Indian Ocean
        end

        @testset "special regions" begin
            # Poles
            @test is_land(π/2, 0) == false                      # North Pole
            @test is_land(-π/2, 0) == true                      # South Pole

            # Deserts
            @test is_land(π/9, deg2rad(5)) == true              # Sahara Desert
            @test is_land(deg2rad(-25), deg2rad(133)) == true   # Australian Outback

            # Around Prime Meridian
            @test is_land(deg2rad(46), deg2rad(2)) == true      # Central France
            @test is_land(0, deg2rad(25)) == true               # Central Africa

            # Around International Date Line
            @test is_land(deg2rad(-41), deg2rad(175)) == true   # New Zealand
            @test is_land(π/4, π) == false                      # Far from land

            # Around Equator
            @test is_land(0, π/9) == true                       # Africa
            @test is_land(0, -2π/3) == false                    # Pacific Ocean
        end

        @testset "return type" begin
            @test is_land((0, 0, 0)) isa Bool  # Tuple input
            @test is_land(0, 0, 0) isa Bool    # Separate arguments
        end

        @testset "consistency" begin
            lat, lon = deg2rad(40.7), deg2rad(-74.0)
            result1 = is_land(lat, lon)
            result2 = is_land(lat, lon)
            @test result1 == result2
        end

        @testset "error handling" begin
            for val in [Inf, -Inf, NaN]
                @test_throws InputError is_land(val, 0)
                @test_throws InputError is_land(0, val)
            end
        end
    end

    @testset "generate_population_center_gses" begin
        @testset "basic generation" begin
            for n in [0, 1, 5, 10, 100]
                gses = generate_population_center_gses(n)
                @test length(gses) == n
                @test count(is_land, gses) >= 0.98 * n
                @test all(in_coord_range, gses)
            end
        end

        @testset "no duplicates" begin
            gses = generate_population_center_gses(20)
            @test length(unique(gses)) == length(gses)
        end

        @testset "min_distance parameter" begin
            gses = generate_population_center_gses(10; min_distance=500)
            @test check_min_distance_constraint(gses, 500)

            gses2 = generate_population_center_gses(10; min_distance=0)
            @test length(gses2) == 10
        end

        @testset "timeout parameter" begin
            gses_small = generate_population_center_gses(10; timeout=5)
            gses_large = generate_population_center_gses(10)    # timeout = Inf
            @test 0 < length(gses_small) <= length(gses_large)
        end

        @testset "verbose parameter" begin
            @test_nowarn generate_population_center_gses(3)

            output = @capture_out generate_population_center_gses(
                5; timeout=3, verbose=true
            )
            @test occursin("Timeout", output)
        end

        @testset "other_gses parameter: builder pattern" begin
            # Returns existing + new stations
            gses = generate_population_center_gses(5; other_gses=[NYC])
            @test length(gses) == 6
            @test is_approx_in(NYC, gses)

            # Multiple existing stations
            existing = [NYC, LONDON, TOKYO, SYDNEY]
            gses = generate_population_center_gses(5; other_gses=existing)
            @test length(gses) == 9
            @test all(gs::GS -> is_approx_in(gs, gses), existing)

            # Empty existing stations edge case
            gses = generate_population_center_gses(5; other_gses=[])
            @test length(gses) == 5

            # Format compatibility
            gses = generate_population_center_gses(5; other_gses=[NYC])
            gses2 = generate_population_center_gses(5; other_gses=[(NYC[1], NYC[2], 0)])
            @test length(gses) == length(gses2) == 6
        end

        @testset "other_gses parameter: distance constraint" begin
            min_dist = 1000
            gses = generate_population_center_gses(
                5; other_gses=[NYC], min_distance=min_dist
            )
            @test check_min_distance_constraint(gses, min_dist)
        end

        @testset "other_gses parameter: timeout constraint" begin
            gses_small = generate_population_center_gses(10; other_gses=[NYC], timeout=5)
            gses_large = generate_population_center_gses(10; other_gses=[NYC], timeout=Inf)
            @test 1 <= length(gses_small) <= length(gses_large)
        end

        @testset "other_gses parameter: consistency" begin
            existing = [NYC, LONDON]
            gses1 = generate_population_center_gses(5; other_gses=existing)
            gses2 = generate_population_center_gses(5; other_gses=existing)
            @test gses1 == gses2
        end

        @testset "other_gses parameter: removes duplicates" begin
            existing = generate_population_center_gses(1)
            gses = generate_population_center_gses(5; other_gses=existing, min_distance=0.1)
            @test length(unique(gses)) == length(gses) == 6
        end

        @testset "return type" begin
            gses = generate_population_center_gses(5)
            @test gses isa Vector
            @test all(gs::GS -> length(gs) == 2, gses)
            @test all(gs::GS -> all(isreal, gs), gses)
        end

        @testset "consistency" begin
            gses1 = generate_population_center_gses(5)
            gses2 = generate_population_center_gses(5)
            @test gses1 == gses2
        end

        @testset "error handling" begin
            for val in [-Inf, -5, NaN, Inf]
                @test_throws InputError generate_population_center_gses(val)
                @test_throws InputError generate_population_center_gses(5; min_distance=val)
            end
            for val in [-Inf, -5, NaN]
                @test_throws InputError generate_population_center_gses(5; timeout=val)
            end
            for val in [-Inf, NaN, Inf]
                for gses in [(val, 0, 0), (0, val, 0), (0, 0, val)]
                    @test_throws InputError generate_population_center_gses(
                        5; other_gses=[gses]
                    )
                end
            end
        end
    end

    @testset "generate_population_center_gses (:rows)" begin
        @testset "basic generation and return type" begin
            for n in [0, 1, 5, 10, 100]
                rows = generate_population_center_gses(Val(:rows), n)
                @test rows isa Vector
                @test length(rows) == n
                if n > 0
                    @test all(r -> length(r) == 11, rows)
                end
            end
        end

        @testset "min_distance parameter" begin
            min_dist = 500
            rows = generate_population_center_gses(Val(:rows), 10; min_distance=min_dist)
            gses = [(deg2rad(parse(Float64, row.lat)), deg2rad(parse(Float64, row.lng))) for row in rows]
            @test check_min_distance_constraint(gses, min_dist)
        end

        @testset "timeout parameter" begin
            rows_small = generate_population_center_gses(Val(:rows), 100; timeout=5)
            rows_large = generate_population_center_gses(Val(:rows), 100)   # timeout = Inf
            @test 0 < length(rows_small) <= length(rows_large)
        end

        @testset "verbose parameter" begin
            @test_nowarn generate_population_center_gses(Val(:rows), 3)

            output = @capture_out generate_population_center_gses(
                Val(:rows), 5; timeout=3, verbose=true
            )
            @test occursin("Timeout", output)
        end

        @testset "consistency" begin
            rows1 = generate_population_center_gses(Val(:rows), 5)
            rows2 = generate_population_center_gses(Val(:rows), 5)
            @test length(rows1) == length(rows2)
        end

        @testset "error handling" begin
            for val in [-Inf, -5, NaN, Inf]
                @test_throws InputError generate_population_center_gses(Val(:rows), val)
                @test_throws InputError generate_population_center_gses(
                    Val(:rows), 5; min_distance=val
                )
            end
            for val in [-Inf, -5, NaN]
                @test_throws InputError generate_population_center_gses(
                    Val(:rows), 5; timeout=val)
            end
        end
    end

    @testset "generate_city_gses" begin
        @testset "basic generation" begin
            gses = generate_city_gses(["New York", "London", "Tokyo", "Sydney"])
            @test length(gses) == 4
            @test all(is_land, gses)
            @test all(in_coord_range, gses)

            gses2 = generate_city_gses(Vector{String}())
            @test length(gses2) == 0
        end

        @testset "unknown city" begin
            gses = generate_city_gses(["New York", "XyzUnknownCity123", "London"])
            @test length(gses) == 2
        end

        @testset "case sensitivity" begin
            gses_lower = generate_city_gses(["new york", "london"])
            gses_upper = generate_city_gses(["NEW YORK", "LONDON"])
            gses_mixed = generate_city_gses(["New York", "London"])
            @test gses_lower == gses_upper == gses_mixed
        end

        @testset "duplicate city names" begin
            city_names = ["Paris", "Paris", "London"]
            gses = generate_city_gses(city_names)
            @test length(gses) == 2
        end

        @testset "special characters" begin
            city_names = ["São Paulo", "Montréal"]
            gses = generate_city_gses(city_names)
            @test length(gses) == 2
        end

        @testset "no duplicates" begin
            city_names = ["New York", "London", "Tokyo"]
            gses = generate_city_gses(city_names)
            @test length(unique(gses)) == length(gses)
        end

        @testset "verbose parameter" begin
            @test_nowarn generate_city_gses(["New York", "UnknownCity123"])

            output = @capture_out generate_city_gses(
                ["New York", "UnknownCity123"]; other_gses=[NYC], verbose=true
            )
            @test occursin("unknowncity123", output)
            @test occursin("Duplicate", output)
        end

        @testset "other_gses parameter: builder pattern" begin
            # Returns existing + new stations
            gses = generate_city_gses(["London", "Tokyo", "Sydney"]; other_gses=[NYC])
            @test length(gses) == 4
            @test is_approx_in(NYC, gses)

            # Multiple existing stations
            gses = generate_city_gses(["London", "Tokyo"]; other_gses=[NYC, SYDNEY])
            @test length(gses) == 4
            @test all(gs::GS -> is_approx_in(gs, gses), [NYC, SYDNEY])

            # Empty existing stations edge case
            gses = generate_city_gses(["London", "Tokyo"]; other_gses=[])
            @test length(gses) == 2

            # Format compatibility
            gses = generate_city_gses(["London", "Tokyo"]; other_gses=[NYC])
            gses2 = generate_city_gses(
                ["London", "Tokyo"]; other_gses=[(NYC[1], NYC[2], 0)]
            )
            @test length(gses) == length(gses2) == 3
            @test all(gs::GS -> is_approx_in(gs, gses2), gses)
        end

        @testset "other_gses parameter: consistency" begin
            existing = [NYC, LONDON]
            gses1 = generate_city_gses(["Tokyo", "Sydney"]; other_gses=existing)
            gses2 = generate_city_gses(["Tokyo", "Sydney"]; other_gses=existing)
            @test gses1 == gses2
        end

        @testset "other_gses parameter: removes duplicates" begin
            existing = generate_city_gses(["New York"])
            gses = generate_city_gses(["New York", "London", "Tokyo"]; other_gses=existing)
            @test length(gses) == 3
        end

        @testset "return type" begin
            gses = generate_city_gses(["New York", "London", "Tokyo"])
            @test gses isa Vector
            @test all(gs::GS -> length(gs) == 2, gses)
            @test all(gs::GS -> all(isreal, gs), gses)
        end

        @testset "consistency" begin
            city_names = ["New York", "London", "Tokyo"]
            gses1 = generate_city_gses(city_names)
            gses2 = generate_city_gses(city_names)
            @test gses1 == gses2
        end

        @testset "error handling" begin
            for val in [-Inf, NaN, Inf]
                for gses in [(val, 0, 0), (0, val, 0), (0, 0, val)]
                    @test_throws InputError generate_city_gses(
                        ["New York"]; other_gses=[gses]
                    )
                end
            end
        end
    end

    @testset "generate_random_gses" begin
        @testset "basic generation" begin
            for n in [0, 1, 5, 10, 100]
                gses = generate_random_gses(n; seed=42)
                @test length(gses) == n
                @test all(in_coord_range, gses)
            end
        end

        @testset "distribution" begin
            gses = generate_random_gses(100; seed=42)
            lats = [lat for (lat, lon) in gses]
            lons = [lon for (lat, lon) in gses]
            @test minimum(lats) < 0 && maximum(lats) > 0
            @test minimum(lons) < 0 && maximum(lons) > 0
        end

        @testset "no duplicates" begin
            gses = generate_random_gses(20)
            @test length(unique(gses)) == length(gses)
        end

        @testset "force_land parameter" begin
            gses = generate_random_gses(20; force_land=true, seed=42)
            @test all(is_land, gses)
        end

        @testset "min_distance parameter" begin
            gses = generate_random_gses(20; min_distance=1000, seed=42)
            @test check_min_distance_constraint(gses, 1000)
        end

        @testset "timeout parameter" begin
            gses = generate_random_gses(100; timeout=10, seed=42)
            @test 1 <= length(gses) <= 10

            gses = generate_random_gses(
                5; timeout=100, force_land=true, min_distance=500, seed=42
            )
            @test length(gses) <= 5
            @test all(is_land, gses)

            gses = generate_random_gses(50; force_land=true, timeout=25, seed=42)
            @test length(gses) <= 25
            @test all(is_land, gses)

            output = @capture_out generate_random_gses(
                50; force_land=true, timeout=25, verbose=true, seed=42
            )
            @test occursin("Timeout", output)
        end

        @testset "seed parameter" begin
            gses1 = generate_random_gses(10; seed=12345)
            gses2 = generate_random_gses(10; seed=12345)
            @test gses1 == gses2
            # @test same_gses(gses1, gses2)

            gses3 = generate_random_gses(10; seed=11111)
            @test gses1 != gses3
            # @test !same_gses(gses1, gses3)

            gses4 = generate_random_gses(10; seed=false)
            gses5 = generate_random_gses(10; seed=false)
            @test gses4 != gses5
            # @test !same_gses(gses4, gses5)
        end

        @testset "verbose parameter" begin
            @test_nowarn generate_random_gses(5)

            output = @capture_out generate_random_gses(5; verbose=true)
            @test occursin("Generating", output)
        end

        @testset "other_gses parameter: builder pattern" begin
            gses = generate_random_gses(5; other_gses=[NYC], seed=42)
            @test length(gses) == 6
            @test is_approx_in(NYC, gses)

            gses = generate_random_gses(5; other_gses=[NYC, LONDON], seed=42)
            @test length(gses) == 7
        end

        @testset "other_gses parameter: distance constraint" begin
            min_dist = 1000
            gses = generate_random_gses(5; other_gses=[NYC], min_distance=min_dist, seed=42)
            @test check_min_distance_constraint(gses, min_dist)
        end

        @testset "other_gses parameter: timeout constraint" begin
            gses_small = generate_random_gses(10; other_gses=[NYC], timeout=5, seed=42)
            gses_large = generate_random_gses(10; other_gses=[NYC], timeout=Inf, seed=42)
            @test 1 <= length(gses_small) <= length(gses_large)
        end

        @testset "other_gses parameter: force_land constraint" begin
            gses = generate_random_gses(
                5; other_gses=[NYC], force_land=true, min_distance=500, seed=42
            )
            @test all(is_land, gses)
        end

        @testset "other_gses parameter: consistency" begin
            existing = [NYC, LONDON]
            gses1 = generate_random_gses(5; other_gses=existing, seed=42)
            gses2 = generate_random_gses(5; other_gses=existing, seed=42)
            @test gses1 == gses2
        end

        @testset "return type" begin
            gses = generate_random_gses(5)
            @test gses isa Vector
            @test all(gs::GS -> length(gs) == 2, gses)
            @test all(gs::GS -> all(isreal, gs), gses)
        end

        @testset "error handling" begin
            for val in [-Inf, -5, NaN, Inf]
                @test_throws InputError generate_random_gses(val)
                @test_throws InputError generate_random_gses(5; min_distance=val)
            end
            for val in [-Inf, -5, NaN]
                @test_throws InputError generate_random_gses(5; timeout=val)
            end
            for val in [-Inf, NaN, Inf]
                for gses in [(val, 0, 0), (0, val, 0), (0, 0, val)]
                    @test_throws InputError generate_random_gses(5; other_gses=[gses])
                end
            end
        end
    end

    @testset "generate_equispaced_gses" begin
        @testset "odd n values" begin
            for n in [1, 3, 5, 11, 21, 101]
                gses = generate_equispaced_gses(n)
                @test length(gses) == n
                @test all(in_coord_range, gses)
            end
        end

        @testset "even n handling" begin
            for n in [2, 4, 6, 8, 10]
                @test_throws ErrorException generate_equispaced_gses(n)

                gses = generate_equispaced_gses(n; fail_on_even=false)
                @test length(gses) == n + 1
            end

            gses = generate_equispaced_gses(7; fail_on_even=false)
            @test length(gses) == 7  # Should not increment odd numbers
        end

        @testset "no duplicates" begin
            gses = generate_equispaced_gses(25)
            @test length(unique(gses)) == length(gses)
        end

        @testset "distribution" begin
            gses = generate_equispaced_gses(25)
            lats = [lat for (lat, lon) in gses]
            lons = [lon for (lat, lon) in gses]
            @test maximum(lats) - minimum(lats) > π/2
            @test minimum(lats) < 0 && maximum(lats) > 0
            @test minimum(lons) < 0 && maximum(lons) > 0

            gses = generate_equispaced_gses(17)
            avg_dist = mean(
                haversine_for_tests(gses[i]..., gses[j]...) for i in eachindex(gses) for
                j in (i+1):length(gses)
            )
            @test avg_dist > 1000  # More than 1000 km for 17 points
        end

        @testset "return type" begin
            gses = generate_equispaced_gses(5)
            @test gses isa Vector
            @test all(gs::GS -> length(gs) == 2, gses)
            @test all(gs::GS -> all(isreal, gs), gses)
        end

        @testset "consistency" begin
            gses1 = generate_equispaced_gses(11)
            gses2 = generate_equispaced_gses(11)
            @test gses1 == gses2
            # @test same_gses(gses1, gses2)
        end

        @testset "error handling" begin
            for val in [-Inf, -5, NaN, Inf]
                @test_throws InputError generate_equispaced_gses(val)
            end
        end
    end

    @testset "generate_grid_gses" begin
        @testset "basic parameters" begin
            gses = generate_grid_gses()
            @test length(gses) > 0

            gses_dense = generate_grid_gses(equatorial_distance_km=200)
            gses_sparse = generate_grid_gses(equatorial_distance_km=800)
            @test length(gses_dense) > length(gses_sparse)
        end

        @testset "coordinate ranges" begin
            gses = generate_grid_gses(equatorial_distance_km=200)
            sample_size = min(50, length(gses))
            sample_indices = sample(1:length(gses), sample_size; replace=false)
            @test all(i -> in_coord_range(gses[i]), sample_indices)

            near_equator = count(gs::GS -> abs(gs[1]) < π/6, gses)
            near_poles = count(gs::GS -> abs(gs[1]) > 5*π/12, gses)
            @test near_equator > near_poles
        end

        @testset "global coverage" begin
            gses = generate_grid_gses(equatorial_distance_km=600)
            lons = [lon for (lat, lon) in gses]
            lats = [lat for (lat, lon) in gses]
            @test minimum(lons) < -π/2 && maximum(lons) > π/2
            @test minimum(lats) < -π/4 && maximum(lats) > π/4
        end

        @testset "no duplicates" begin
            gses = generate_grid_gses()
            @test length(unique(gses)) == length(gses)
        end

        @testset "force_land parameter" begin
            gses_all = generate_grid_gses(equatorial_distance_km=500)
            gses_land = generate_grid_gses(equatorial_distance_km=500; force_land=true)
            @test length(gses_land) <= length(gses_all)

            sample_size = min(50, length(gses_land))
            sample_indices = sample(1:length(gses_land), sample_size; replace=false)
            @test all(i -> is_land(gses_land[i]), sample_indices)
        end

        @testset "discount_func parameter" begin
            discount_identity = lat_rad -> 0
            gses = generate_grid_gses(
                equatorial_distance_km=500; discount_func=discount_identity
            )
            @test length(gses) > 0

            discount_poles = lat_rad -> abs(lat_rad) / (π/2)
            gses = generate_grid_gses(
                equatorial_distance_km=500; discount_func=discount_poles
            )

            near_poles = count(gs::GS -> abs(gs[1]) > 5*π/12, gses)
            near_equator = count(gs::GS -> abs(gs[1]) < π/8, gses)
            @test near_equator >= near_poles
        end

        @testset "equatorial_distance_km parameter" begin
            gses = generate_grid_gses(equatorial_distance_km=500)
            near_equator = count(gs::GS -> abs(gs[1]) < 0.1, gses)
            @test near_equator > 0

            gses = generate_grid_gses(equatorial_distance_km=2000)
            @test length(gses) >= 1

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

        @testset "other_gses parameter: builder pattern" begin
            existing = [NYC, LONDON]
            gses = generate_grid_gses(equatorial_distance_km=500; other_gses=existing)
            @test length(gses) >= length(existing)
            @test all(gs::GS -> is_approx_in(gs, gses), existing)
        end

        @testset "other_gses parameter: with force_land" begin
            gses = generate_grid_gses(
                equatorial_distance_km=500; other_gses=[NYC], force_land=true
            )
            @test all(is_land, gses)
        end

        @testset "other_gses parameter: consistency" begin
            gses1 = generate_grid_gses(equatorial_distance_km=500; other_gses=[NYC])
            gses2 = generate_grid_gses(equatorial_distance_km=500; other_gses=[NYC])
            @test gses1 == gses2
        end

        @testset "return type" begin
            gses = generate_grid_gses()
            sample_size = min(20, length(gses))
            sample_indices = sample(1:length(gses), sample_size; replace=false)

            @test all(i -> gses[i] isa Tuple && length(gses[i]) == 2, sample_indices)
            @test all(i -> all(isreal, gses[i]), sample_indices)
        end

        @testset "consistency" begin
            gses1 = generate_grid_gses(equatorial_distance_km=400)
            gses2 = generate_grid_gses(equatorial_distance_km=400)
            @test gses1 == gses2
            # @test same_gses(gses1, gses2)
        end

        @testset "error handling" begin
            for val in [-Inf, -5, NaN, Inf]
                @test_throws InputError generate_grid_gses(equatorial_distance_km=val)
            end
        end
    end
end