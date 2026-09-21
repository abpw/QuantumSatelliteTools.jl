using Dates: DateTime
using Test: @testset, @test, @test_throws, @test_nowarn, @test_broken

using StaticArrays: SVector
using SatelliteToolboxBase: JD_J2000, Ellipsoid, OrbitStateVector, date_to_jd
using SatelliteToolboxTle: read_tle
using SatelliteToolboxPropagators: Propagators
using SatelliteAnalysis: sun_position_mod

using QuantumSatelliteTools
using QuantumSatelliteTools.constants
using QuantumSatelliteTools.types
using QuantumSatelliteTools.helpers
using QuantumSatelliteTools.AstronomyGeometry

include("test_tles.jl")

@testset "AstronomyGeometry.jl" begin
    @testset "norm" begin
        @testset "zero vector" begin
            @test norm((0, 0, 0)) == 0
        end

        @testset "unit vectors" begin
            @test isapprox(norm((1, 0, 0)), 1; atol=1e-9)
            @test isapprox(norm((0, 1, 0)), 1; atol=1e-9)
            @test isapprox(norm((0, 0, 1)), 1; atol=1e-9)
        end

        @testset "3-4-5 triangle" begin
            @test isapprox(norm((3, 4, 0)), 5; atol=1e-9)
            @test isapprox(norm((-3, -4, 0)), 5; atol=1e-9)
        end

        @testset "diagonal of unit cube" begin
            @test isapprox(norm((1, 1, 1)), sqrt(3); atol=1e-9)
        end

        @testset "very small values" begin
            @test isapprox(norm((1e-9, 1e-9, 1e-9)), sqrt(3) * 1e-9; atol=1e-9)
        end

        @testset "very large values" begin
            @test isapprox(norm((1e10, 1e10, 1e10)), sqrt(3) * 1e10; rtol=1e-9)
        end

        @testset "Tuple input" begin
            @test isapprox(norm((3, 4, 0)), 5; atol=1e-9)
        end

        @testset "SVector input" begin
            @test isapprox(norm(SVector(3, 4, 0)), 5; atol=1e-9)
        end

        @testset "wrong dimension" begin
            @test_throws MethodError norm((1, 2))
            @test_throws MethodError norm((1, 2, 3, 4))
        end

        @testset "wrong type" begin
            @test_throws MethodError norm([1, 2, 3])
            @test_throws MethodError norm(("a", "b", "c"))
            @test_throws MethodError norm((1 + 1im, 2, 3))
        end
    end

    # Testing conversions in case implementation changes, currently external
    @testset "geodetic_to_ecef" begin
        @testset "origin" begin
            ecef = geodetic_to_ecef((0, 0, 0))
            expected_ecef = [SEMIMAJOR_RADIUS, 0, 0]

            @test ecef isa AbstractVector
            @test all(isreal.(ecef))
            @test all(isfinite.(ecef))
            @test length(ecef) == 3

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "equator, 90° east" begin
            ecef = geodetic_to_ecef((0, π/2, 0))
            expected_ecef = [0, SEMIMAJOR_RADIUS, 0]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "equator, 90° west" begin
            ecef = geodetic_to_ecef((0, -π/2, 0))
            expected_ecef = [0, -SEMIMAJOR_RADIUS, 0]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "equator, 180°" begin
            ecef = geodetic_to_ecef((0, π, 0))
            expected_ecef = [-SEMIMAJOR_RADIUS, 0, 0]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "north pole" begin
            ecef = geodetic_to_ecef((π/2, 0, 0))
            expected_ecef = [0, 0, SEMIMINOR_RADIUS]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "south pole" begin
            ecef = geodetic_to_ecef((-π/2, 0, 0))
            expected_ecef = [0, 0, -SEMIMINOR_RADIUS]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "near north pole" begin
            ecef = geodetic_to_ecef((deg2rad(89.9999999), 0, 0))
            expected_ecef = [0.011169397956, 0, 6_356_752.314245]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "North pole with nonzero longitude" begin
            ecef = geodetic_to_ecef((π/2, deg2rad(137), 0))
            expected_ecef = [0, 0, SEMIMINOR_RADIUS]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "shallow negative altitude" begin
            ecef = geodetic_to_ecef((0, 0, -1_000))
            expected_ecef = [SEMIMAJOR_RADIUS - 1_000, 0, 0]

            @test(all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9)))
        end

        @testset "negative semimajor altitude" begin
            ecef = geodetic_to_ecef((0, 0, -SEMIMAJOR_RADIUS))
            expected_ecef = [0, 0, 0]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9))
        end

        @testset "GEO altitude" begin
            ecef = geodetic_to_ecef((0, 0, 35_786_000))
            expected_ecef = [42_164_137, 0, 0]

            @test all(isapprox.(ecef, expected_ecef; atol=1e-9, rtol=1e-9))
        end

        @testset "2-tuple input (no height)" begin
            ecef = geodetic_to_ecef((
                deg2rad(42.394205296564806), deg2rad(-72.52761607116234)
            ))
            expected_ecef = [1_416_435.823714, -4_499_938.662610, 4_278_031.082705]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-5))
        end

        @testset "northeast" begin
            ecef = geodetic_to_ecef((π/4, π/4, 0))
            expected_ecef = [3_194_419.145061, 3_194_419.145061, 4_487_348.408866]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-9))
        end

        @testset "northeast, 1000 m height" begin
            ecef = geodetic_to_ecef((π/4, π/4, 1000.0))
            expected_ecef = [3_194_919.145061, 3_194_919.145061, 4_488_055.515647]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-9))
        end

        @testset "southwest" begin
            ecef = geodetic_to_ecef((-π/4, -π/4, 0.0))
            expected_ecef = [3_194_419.145061, -3_194_419.145061, -4_487_348.408866]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-9))
        end

        @testset "DC, 100 m height" begin
            ecef = geodetic_to_ecef((deg2rad(39), deg2rad(-77), 100.0))
            expected_ecef = [1_116_523.199853, -4_836_193.303241, 3_992_379.954791]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-9))
        end

        @testset "180° = -180°" begin
            ecef_180 = geodetic_to_ecef((0.0, π, 0.0))
            ecef_neg_180 = geodetic_to_ecef((0.0, -π, 0.0))

            @test all(isapprox.(ecef_180, ecef_neg_180; atol=1e-8, rtol=1e-9))
        end

        @testset "359° = -1°" begin
            ecef_359 = geodetic_to_ecef((0.0, deg2rad(359), 0.0))
            ecef_neg_1 = geodetic_to_ecef((0.0, deg2rad(-1), 0.0))

            @test all(isapprox.(ecef_359, ecef_neg_1; rtol=1e-9))
        end

        @testset "450° = 90°" begin
            ecef_450 = geodetic_to_ecef((0.0, 5π/2, 0.0))
            ecef_90 = geodetic_to_ecef((0.0, π/2, 0.0))

            @test all(isapprox.(ecef_450, ecef_90; atol=1e-8, rtol=1e-9))
        end

        @testset "infinite values" begin
            @test_throws ArgumentError geodetic_to_ecef((Inf, 0.0, 0.0))
            @test_throws ArgumentError geodetic_to_ecef((0.0, Inf, 0.0))
            @test_throws ArgumentError geodetic_to_ecef((0.0, 0.0, Inf))
        end

        @testset "NaN values" begin
            @test_throws ArgumentError geodetic_to_ecef((NaN, 0.0, 0.0))
            @test_throws ArgumentError geodetic_to_ecef((0.0, NaN, 0.0))
            @test_throws ArgumentError geodetic_to_ecef((0.0, 0.0, NaN))
        end

        @testset "different ellipsoid" begin
            ecef = geodetic_to_ecef(
                (deg2rad(39), deg2rad(-77), 100.0),
                ellipsoid=Ellipsoid(6378137, 1 / 298.257222101),
            )
            expected_ecef = [1_116_523.199860, -4_836_193.303272, 3_992_379.954685]

            @test all(isapprox.(ecef, expected_ecef; rtol=1e-9))
        end
    end

    @testset "ecef_to_geodetic" begin
        @testset "geodetic origin" begin
            geodetic = ecef_to_geodetic((SEMIMAJOR_RADIUS, 0.0, 0.0))
            expected_geodetic = [0, 0, 0]

            @test geodetic isa Tuple
            @test all(isreal.(geodetic))
            @test all(isfinite.(geodetic))
            @test length(geodetic) == 3

            @test all(isapprox.([geodetic...], expected_geodetic; atol=1e-9, rtol=1e-9))
        end

        @testset "equator, 90° east" begin
            geodetic = [ecef_to_geodetic((0.0, SEMIMAJOR_RADIUS, 0.0))...]
            expected_geodetic = [0, π/2, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-9))
        end

        @testset "equator, 90° west" begin
            geodetic = [ecef_to_geodetic((0.0, -SEMIMAJOR_RADIUS, 0.0))...]
            expected_geodetic = [0, -π/2, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-9))
        end

        @testset "equator, 180°" begin
            geodetic = [ecef_to_geodetic((-SEMIMAJOR_RADIUS, 0.0, 0.0))...]
            expected_geodetic = [0, π, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-9))
        end

        @testset "north pole" begin
            geodetic = [ecef_to_geodetic((0.0, 0.0, SEMIMINOR_RADIUS))...]
            expected_geodetic = [π/2, 0, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "south pole" begin
            geodetic = [ecef_to_geodetic((0.0, 0.0, -SEMIMINOR_RADIUS))...]
            expected_geodetic = [-π/2, 0, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "near north pole" begin
            geodetic = [ecef_to_geodetic((0.009673, 0.005585, 6_356_752.314245))...]
            expected_geodetic = [deg2rad(89.9999999), deg2rad(30.00129205), 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-5))
        end

        @testset "shallow negative altitude" begin
            geodetic = [ecef_to_geodetic((SEMIMAJOR_RADIUS - 1_000, 0.0, 0.0))...]
            expected_geodetic = [0, 0, -1_000]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4, rtol=1e-4))
        end

        @testset "geocenter" begin
            geodetic1 = [ecef_to_geodetic((0.0, 0.0, 0.0))...]
            geodetic2 = [ecef_to_geodetic((0.0, 0.0, 0.0))...]

            @test_nowarn ecef_to_geodetic((0.0, 0.0, 0.0))
            @test all(isfinite.(geodetic1))
            @test all(isapprox.(geodetic1, geodetic2; atol=1e-4))
        end

        @testset "GEO altitude" begin
            geodetic = [ecef_to_geodetic((42_164_137, 0.0, 0.0))...]
            expected_geodetic = [0, 0, 35_786_000]

            @test all(isapprox.(geodetic, expected_geodetic; rtol=1e-4))
        end

        @testset "northeast" begin
            geodetic = [
                ecef_to_geodetic((3_194_419.145061, 3_194_419.145061, 4_487_348.408866))...
            ]
            expected_geodetic = [π/4, π/4, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "northeast, 1000 m height" begin
            geodetic = [
                ecef_to_geodetic((3_194_919.145061, 3_194_919.145061, 4_488_055.515647))...
            ]
            expected_geodetic = [π/4, π/4, 1000]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "southwest" begin
            geodetic = [
                ecef_to_geodetic((
                    3_194_419.145061, -3_194_419.145061, -4_487_348.408866
                ))...,
            ]
            expected_geodetic = [-π/4, -π/4, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "DC, 100 m height" begin
            geodetic = [
                ecef_to_geodetic((1_116_523.199853, -4_836_193.303241, 3_992_379.954791))...
            ]
            expected_geodetic = [deg2rad(39), deg2rad(-77), 100]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "positive signed zero" begin
            geodetic = [ecef_to_geodetic((-SEMIMAJOR_RADIUS, +0.0, 0.0))...]
            expected_geodetic = [0, π, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "negative signed zero" begin
            geodetic = [ecef_to_geodetic((-SEMIMAJOR_RADIUS, -0.0, 0.0))...]
            expected_geodetic = [0, -π, 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "small positive perturbation around zero" begin
            geodetic = [ecef_to_geodetic((-SEMIMAJOR_RADIUS, 0.001, 0.0))...]
            expected_geodetic = [0, deg2rad(179.9999999910), 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "small negative perturbation around zero" begin
            geodetic = [ecef_to_geodetic((-SEMIMAJOR_RADIUS, -0.001, 0.0))...]
            expected_geodetic = [0, deg2rad(-179.9999999910), 0]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "infinite values" begin
            @test_throws ArgumentError ecef_to_geodetic((Inf, 0.0, 0.0))
            @test_throws ArgumentError ecef_to_geodetic((0.0, Inf, 0.0))
            @test_throws ArgumentError ecef_to_geodetic((0.0, 0.0, Inf))
        end

        @testset "NaN values" begin
            @test_throws ArgumentError ecef_to_geodetic((NaN, 0.0, 0.0))
            @test_throws ArgumentError ecef_to_geodetic((0.0, NaN, 0.0))
            @test_throws ArgumentError ecef_to_geodetic((0.0, 0.0, NaN))
        end

        @testset "different ellipsoid" begin
            geodetic = [
                ecef_to_geodetic(
                    (1_116_523.199860, -4_836_193.303272, 3_992_379.954685),
                    ellipsoid=Ellipsoid(6378137, 1 / 298.257222101),
                )...,
            ]
            expected_geodetic = [deg2rad(39), deg2rad(-77), 100]

            @test all(isapprox.(geodetic, expected_geodetic; atol=1e-4))
        end

        @testset "longitude range convention" begin
            points = [
                (SEMIMAJOR_RADIUS, 0.0, 0.0),
                (0.0, SEMIMAJOR_RADIUS, 0.0),
                (0.0, -SEMIMAJOR_RADIUS, 0.0),
                (-SEMIMAJOR_RADIUS, 0.0, 0.0),
                # (-SEMIMAJOR_RADIUS, -0.0, 0.0), - broken
            ]
            longitudes = [ecef_to_geodetic(point)[2] for point in points]

            @test all(-π < lon <= π for lon in longitudes)
        end
    end

    @testset "geodetic_to_ecef & ecef_to_geodetic roundtrip" begin
        test_points = [
            (0.0, 0.0, 0.0),
            (π/2, 0.0, 0.0),
            (-π/2, 0.0, 0.0),
            (π/4, π/3, 1000.0),
            (π/4, -2π/3, 500_000.0),
            (deg2rad(42.39466434352798), deg2rad(-72.5270869156386), 80.0),
        ]
        roundtrips = [ecef_to_geodetic(geodetic_to_ecef(point)) for point in test_points]

        @test all(
            all(isapprox.(point, roundtrip_point; rtol=1e-4)) for
            (point, roundtrip_point) in zip(test_points, roundtrips)
        )
    end

    @testset "eci_to_ecef" begin
        @testset "OrbitStateVector input" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            ecef = eci_to_ecef(sv)

            @test length(ecef) == 3
            @test all(isreal.(ecef))
            @test all(isfinite.(ecef))
            @test ecef isa AbstractVector

            r = √(sum(ecef .^ 2))

            @test r > SEMIMAJOR_RADIUS
            @test r < SEMIMAJOR_RADIUS + 1_000_000
        end

        @testset "Propagator input with time" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            ecef = eci_to_ecef(sat; time=sat.sgp4d.epoch)

            @test length(ecef) == 3
            @test all(isreal.(ecef))
            @test all(isfinite.(ecef))
            @test ecef isa AbstractVector

            r = √(sum(ecef .^ 2))

            @test r > SEMIMAJOR_RADIUS
        end

        @testset "consistency across input formats" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            ecef_sv = eci_to_ecef(sv; time=sat.sgp4d.epoch)
            ecef_tuple = eci_to_ecef(Tuple(sv.r); time=sat.sgp4d.epoch)
            ecef_xyz = eci_to_ecef(sv.r[1], sv.r[2], sv.r[3]; time=sat.sgp4d.epoch)

            @test all(isapprox.(ecef_sv, ecef_tuple; rtol=1e-10))
            @test all(isapprox.(ecef_sv, ecef_xyz; rtol=1e-10))
        end

        @testset "handles NaN/Inf gracefully" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            @test_throws ArgumentError eci_to_ecef(
                sv.r[1], Inf, sv.r[3]; time=sat.sgp4d.epoch
            )
        end

        @testset "different times produce different results" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            ecef1 = eci_to_ecef(Tuple(sv.r); time=sat.sgp4d.epoch)
            ecef2 = eci_to_ecef(Tuple(sv.r); time=sat.sgp4d.epoch + 3600)  # 1 hour later

            @test !all(isapprox.(ecef1, ecef2; rtol=1e-6))
        end
    end

    @testset "eci_to_geodetic" begin
        @testset "OrbitStateVector input" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            geodetic = eci_to_geodetic(sv)

            @test all(isfinite.(geodetic))
            @test all(isreal.(geodetic))
            @test length(geodetic) == 3
            @test geodetic isa Tuple

            @test -π/2 <= geodetic[1] <= π/2
            @test -π <= geodetic[2] <= π
            @test geodetic[3] > 400_000    # ISS altitude ~400+ km
        end

        @testset "Propagator input with time" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            lat, lon, h = eci_to_geodetic(sat; time=sat.sgp4d.epoch)

            @test -π/2 <= lat <= π/2
            @test -π <= lon <= π
            @test h > 400_000    # ISS altitude ~400+ km
        end

        @testset "consistency across input formats" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            geodetic_sv = [eci_to_geodetic(sv)...]
            geodetic_tuple = [eci_to_geodetic(Tuple(sv.r); time=sat.sgp4d.epoch)...]
            geodetic_xyz = [
                eci_to_geodetic(sv.r[1], sv.r[2], sv.r[3]; time=sat.sgp4d.epoch)...
            ]

            @test all(isapprox.(geodetic_sv, geodetic_tuple; rtol=1e-6))
            @test all(isapprox.(geodetic_sv, geodetic_xyz; rtol=1e-6))
        end

        @testset "infinite values" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            @test_throws ArgumentError eci_to_geodetic(
                Inf, sv.r[2], sv.r[3]; time=sat.sgp4d.epoch
            )
        end

        @testset "NaN values" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)

            @test_throws ArgumentError eci_to_geodetic(
                NaN, sv.r[2], sv.r[3]; time=sat.sgp4d.epoch
            )
        end

        @testset "validates altitude is positive" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            h = eci_to_geodetic(sv)[3]

            @test h >= 0
        end

        @testset "validates latitude range" begin
            tle = read_tle(TEST_TLES[:ISS]; verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            lat = eci_to_geodetic(sv)[1]

            @test lat >= -π/2 && lat <= π/2
        end
    end

    @testset "ellipsoid_line_intersection" begin
        @testset "line through center" begin
            p1, p2 = (2, 0, 0), (-2, 0, 0)
            xs = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs) == 2
            @test all(isapprox(sqrt(sum(Tuple(x) .^ 2)), 1; atol=1e-9) for x in xs)
        end

        @testset "line misses ellipsoid" begin
            p1, p2 = (10, 10, 0), (10, 20, 0)
            xs = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs) == 0
        end

        @testset "tangent line" begin
            p1, p2 = (1, 5, 0), (1, -5, 0)
            xs = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs) == 1
        end

        @testset "oblate ellipsoid" begin
            xs = ellipsoid_line_intersection(
                SEMIMAJOR_RADIUS,
                SEMIMINOR_RADIUS,
                (SEMIMAJOR_RADIUS * 2, 0, 0),
                (-SEMIMAJOR_RADIUS * 2, 0, 0),
            )

            @test length(xs) == 2
            @test all(isapprox(x[3], 0; atol=1e-9) for x in xs)
        end

        @testset "no points between" begin
            p1, p2 = (2, 0, 0), (3, 0, 0)
            xs_between = ellipsoid_line_intersection(1, 1, p1, p2; only_between=true)
            xs_all = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs_between) == 0
            @test length(xs_all) == 2
        end

        @testset "one point between" begin
            p1, p2 = (2, 0, 0), (0.5, 0, 0)
            xs_between = ellipsoid_line_intersection(1, 1, p1, p2; only_between=true)
            xs_all = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs_between) == 1
            @test length(xs_all) == 2
        end

        @testset "two points between" begin
            p1, p2 = (-2, 0, 0), (2, 0, 0)
            xs_between = ellipsoid_line_intersection(1, 1, p1, p2; only_between=true)
            xs_all = ellipsoid_line_intersection(1, 1, p1, p2)

            @test length(xs_between) == 2
            @test length(xs_all) == 2
        end

        @testset "intersection properties" begin
            major, minor = 2, 1.5
            p1, p2 = (4, 0, 0), (-4, 0, 0)
            xs = ellipsoid_line_intersection(major, minor, p1, p2)

            for x in xs
                @test isapprox(x[2], 0; atol=1e-9)
                @test isapprox(x[3], 0; atol=1e-9)

                ellipsoid_eq = (x[1]^2 / major^2) + (x[2]^2 / minor^2) + (x[3]^2 / minor^2)
                @test isapprox(ellipsoid_eq, 1; atol=1e-9)
            end
        end

        @testset "return type" begin
            p1, p2 = (2, 0, 0), (-2, 0, 0)
            xs = ellipsoid_line_intersection(1, 1, p1, p2)

            @test xs isa Vector
            @test all(x isa Tuple for x in xs)
            @test all(length(x) == 3 for x in xs)
        end

        @testset "number of intersections" begin
            test_cases = [
                (1, 1, (2, 0, 0), (-2, 0, 0)),
                (1, 1, (10, 10, 0), (10, 20, 0)),
                (
                    SEMIMAJOR_RADIUS,
                    SEMIMINOR_RADIUS,
                    (SEMIMAJOR_RADIUS*2, 0.0, 0.0),
                    (-SEMIMAJOR_RADIUS*2, 0.0, 0.0),
                ),  # Earth
            ]
            xss = [ellipsoid_line_intersection(args...) for args in test_cases]

            @test all(length(xs) in [0, 1, 2] for xs in xss)
        end
    end

    @testset "gs_gs_distance" begin
        @testset "same location" begin
            gs, gs_mid = (0, 0, 0), (π/4, -2π/3, 100.0)
            dist = gs_gs_distance(gs, gs)
            dist_mid = gs_gs_distance(gs_mid, gs_mid)

            @test isapprox(dist, 0; atol=1)
            @test isapprox(dist_mid, 0; atol=1)
        end

        @testset "symmetry" begin
            gs1, gs2 = (0, 0, 0), (π/4, π/2, 0.0)
            dist1 = gs_gs_distance(gs1, gs2)
            dist2 = gs_gs_distance(gs2, gs1)

            @test isapprox(dist1, dist2; rtol=1e-9)
        end

        @testset "non-negative" begin
            test_pairs = [
                ((0, 0, 0), (π/4, π/4, 0.0)),
                ((π/2, 0.0, 0.0), (-π/2, 0.0, 0.0)),
                ((π/4, -2π/3, 100.0), (deg2rad(35), deg2rad(-110), 500.0)),
            ]
            dists = [gs_gs_distance(gs1, gs2) for (gs1, gs2) in test_pairs]

            @test all(dist >= 0 for dist in dists)
        end

        @testset "equator quarter circle" begin
            gs1, gs2 = (0, 0, 0), (0.0, π/2, 0.0)
            dist = gs_gs_distance(gs1, gs2)
            expected_dist = (π/2) * SEMIMAJOR_RADIUS

            @test isapprox(dist, expected_dist; rtol=0.01)
        end

        @testset "explicit format tags" begin
            gs1, gs2 = (0, 0, 0), (π/4, π/2, 0.0)
            dist = gs_gs_distance(gs1, gs2)
            dist_latlon = gs_gs_distance(Val(:lat_lon), gs1, gs2)

            @test isapprox(dist, dist_latlon; rtol=1e-9)
        end

        @testset "ECEF format" begin
            ecef1, ecef2 = (SEMIMAJOR_RADIUS, 0.0, 0.0), (0.0, SEMIMAJOR_RADIUS, 0.0)
            dist = gs_gs_distance(Val(:ECEF), ecef1, ecef2)
            expected_dist = (π/2) * SEMIMAJOR_RADIUS

            @test isapprox(dist, expected_dist; rtol=0.01)
        end

        @testset "with altitude" begin
            gs_sea, gs_alt = (0, 0, 0), (0, 0, 1_000)
            dist = gs_gs_distance(gs_sea, gs_alt)

            @test isapprox(dist, 1_000; rtol=1e-9)
        end

        @testset "two altitude points" begin
            gs1, gs2 = (0, 0, 5_000), (0, 0, 10_000)
            dist = gs_gs_distance(gs1, gs2)

            @test dist > 0
            @test isapprox(dist, 5_000; rtol=1e-3)
        end

        @testset "poles" begin
            gs1, gs2 = (π/2, 0, 0), (-π/2, 0, 0)
            dist = gs_gs_distance(gs1, gs2)
            expected_dist = π * SEMIMINOR_RADIUS

            @test isapprox(dist, expected_dist; rtol=0.01)
        end

        @testset "equator antipode" begin
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (0.0, π, 0.0)

            dist = gs_gs_distance(gs1, gs2)
            expected_dist = π * SEMIMAJOR_RADIUS
            @test isapprox(dist, expected_dist; rtol=0.01)
        end

        @testset "return type" begin
            gs1, gs2 = (0, 0, 0), (π/4, π/4, 0.0)
            dist = gs_gs_distance(gs1, gs2)

            @test dist isa Number
            @test dist >= 0
        end

        @testset "short distances" begin
            gs1, gs2 = (0, 0, 0), (1e-6, 1e-6, 0.0)
            dist = gs_gs_distance(gs1, gs2)

            @test dist > 0
            @test dist < 1000
        end

        @testset "2-tuple input" begin
            gs1, gs2 = (0, 0), (0.0, π/2)
            dist = gs_gs_distance(gs1, gs2)
            expected_dist = (π/2) * SEMIMAJOR_RADIUS

            @test isapprox(dist, expected_dist; rtol=0.01)
        end
    end

    @testset "sat_sat_distance" begin
        @testset "same satellite" begin
            tle = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)
            dist = sat_sat_distance(sat, sat; time=sat.sgp4d.epoch)

            @test dist === nothing || isapprox(dist, 0; atol=1e-9)
        end

        @testset "non-negative when not obstructed" begin
            tle1 = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            tle2 = read_tle(TEST_TLES[:Hubble], verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            dist = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)

            if dist !== nothing
                @test dist >= 0
                @test dist > 1000
            end
        end

        @testset "time parameter" begin
            tle1 = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            tle2 = read_tle(TEST_TLES[:Hubble], verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            dist1 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)
            dist2 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch + 100)

            if dist1 !== nothing && dist2 !== nothing
                @test isa(dist1, Number)
                @test isa(dist2, Number)
            end
        end

        @testset "return type" begin
            tle = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)

            dist = sat_sat_distance(sat, sat; time=sat.sgp4d.epoch)
            @test dist === nothing || dist isa Number
        end

        @testset "LEO to LEO" begin
            tle1 = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            tle2 = read_tle(TEST_TLES[:NOAA18], verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            dist = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)
            if dist !== nothing
                @test dist > 1_000_000 # ISS orbit ~400km, NOAA18 ~800km, min distance ~1000km
                @test dist < 20_000_000
            end
        end

        @testset "symmetry" begin
            tle1 = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            tle2 = read_tle(TEST_TLES[:Hubble], verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            dist12 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)
            dist21 = sat_sat_distance(sat2, sat1; time=sat1.sgp4d.epoch)

            if dist12 !== nothing && dist21 !== nothing
                @test isapprox(dist12, dist21; rtol=1e-9)
            end
        end

        @testset "time dependency over 24 hours" begin
            tle1 = read_tle(TEST_TLES[:ISS], verify_checksum=false)
            tle2 = read_tle(TEST_TLES[:NOAA18], verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            dists = []
            for hours in 0:6:24
                dist = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch + hours/24)
                dist !== nothing && push!(dists, dist)
            end

            if length(dists) > 1
                @test maximum(dists) != minimum(dists)
            end
        end
    end

    @testset "sun_position" begin
        @testset "return type and magnitude" begin
            sun = sun_position(JD_J2000)
            r = sqrt(sum(sun .^ 2))

            @test sun isa AbstractVector
            @test length(sun) == 3
            @test isapprox(r, AU_M; atol=2.5e9)
        end

        @testset "DateTime input" begin
            dt = DateTime(2000, 1, 1, 12, 0, 0)
            sun_dt = sun_position(dt)
            r = sqrt(sum(sun_dt .^ 2))

            @test sun_dt isa AbstractVector
            @test length(sun_dt) == 3
            @test isapprox(r, AU_M; atol=2.5e9)
        end

        @testset "different epochs" begin
            times = [JD_J2000, JD_J2000 + 1, JD_J2000 + 365, JD_J2000 + 365*10]

            for jd in times
                sun = sun_position(jd)
                r = sqrt(sum(sun .^ 2))

                @test isapprox(r, AU_M; atol=2.5e9)
            end
        end

        @testset "orbital motion" begin
            sun1 = sun_position(JD_J2000)
            sun2 = sun_position(JD_J2000 + 90)  # 90 days later

            cos_angle =
                (sun1[1]*sun2[1] + sun1[2]*sun2[2] + sun1[3]*sun2[3]) /
                    (sqrt(sum(sun1 .^ 2)) * sqrt(sum(sun2 .^ 2)))
            angle = acos(clamp(cos_angle, -1, 1))

            @test isapprox(angle, π/2; rtol=0.01)
        end

        @testset "consistency" begin
            sun1 = sun_position(JD_J2000)
            sun2 = sun_position(JD_J2000)

            @test all(isapprox.(sun1, sun2; rtol=1e-9))
        end
    end

    @testset "light_at_point" begin
        @testset "return type" begin
            point = (1e11, 0.0, 0.0)
            result = light_at_point(point, JD_J2000)

            @test result isa LightCondition
            @test result in [sunlight, penumbra, umbra]
        end

        @testset "consistency" begin
            point = (1e11, 0.0, 0.0)

            condition1 = light_at_point(point, JD_J2000)
            condition2 = light_at_point(point, JD_J2000)

            @test condition1 == condition2
        end

        @testset "DateTime input" begin
            point = (1e11, 0.0, 0.0)
            dt = DateTime(2000, 1, 1, 12, 0, 0)

            condition_dt = light_at_point(point, dt)
            condition_jd = light_at_point(point, JD_J2000)

            @test condition_dt in [sunlight, penumbra, umbra]
            @test condition_dt == condition_jd
        end

        @testset "different times give valid results" begin
            point = (1e11, 0.0, 0.0)

            condition1 = light_at_point(point, JD_J2000)
            condition2 = light_at_point(point, JD_J2000 + 100)

            @test condition1 in [sunlight, penumbra, umbra]
            @test condition2 in [sunlight, penumbra, umbra]
        end

        @testset "random points all valid" begin
            for _ in 1:20
                x = randn() * 1e11
                y = randn() * 1e11
                z = randn() * 1e11
                cond = light_at_point((x, y, z), JD_J2000)

                @test cond in [sunlight, penumbra, umbra]
            end
        end

        @testset "specific known cases" begin
            # borrowed from SatelliteAnalysis:light_condition 
            # tested here in case implementation changes
            time = DateTime(2021, 12, 8, 0, 0, 25)

            @testset "sunlight" begin
                point = (1.06746e6, 3.22132e6, -6.2841e6)
                condition = light_at_point(point, time)

                @test condition == sunlight
            end

            @testset "penumbra" begin
                point = (2.52449e6, 4.89401e6, -4.54648e6)
                condition = light_at_point(point, time)

                @test condition == penumbra
            end

            @testset "umbra" begin
                point = (3.9575e6, 5.91988e6, -489577.0)
                condition = light_at_point(point, time)

                @test condition == umbra
            end

            @testset "direct sunlight" begin
                sun_pos = sun_position_mod(date_to_jd(time))
                point = sun_pos / norm(sun_pos) * 7000e3  # 7000 km in direction of Sun
                condition = light_at_point(point, time)

                @test condition == sunlight
            end
        end
    end
end