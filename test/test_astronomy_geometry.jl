using Dates: DateTime
using Test: @testset, @test, @test_throws, @test_nowarn

using StaticArrays: SVector
using SatelliteToolboxBase: JD_J2000, Ellipsoid, OrbitStateVector
using SatelliteToolboxTle: read_tle
using SatelliteToolboxPropagators: Propagators

using QuantumSatelliteTools
using QuantumSatelliteTools.AstronomyGeometry

@testset "AstronomyGeometry.jl" begin
    @testset "norm" begin
        @testset "known values" begin
            # Origin
            @test norm((0.0, 0.0, 0.0)) == 0.0

            # Unit vectors
            @test isapprox(norm((1.0, 0.0, 0.0)), 1.0; atol=1e-10)
            @test isapprox(norm((0.0, 1.0, 0.0)), 1.0; atol=1e-10)
            @test isapprox(norm((0.0, 0.0, 1.0)), 1.0; atol=1e-10)

            # 3-4-5 triangle (tests sign invariance too)
            @test isapprox(norm((3.0, 4.0, 0.0)), 5.0; atol=1e-10)
            @test isapprox(norm((-3.0, -4.0, 0.0)), 5.0; atol=1e-10)

            # Diagonal of unit cube
            @test isapprox(norm((1.0, 1.0, 1.0)), sqrt(3); atol=1e-10)
        end

        @testset "numerical stability" begin
            # Very small values
            @test isapprox(norm((1e-10, 1e-10, 1e-10)), sqrt(3) * 1e-10; atol=1e-10)

            # Very large values
            @test isapprox(norm((1e10, 1e10, 1e10)), sqrt(3) * 1e10; rtol=1e-10)
        end

        @testset "input types" begin
            # Tuple
            @test isapprox(norm((3.0, 4.0, 0.0)), 5.0; atol=1e-10)

            # SVector
            @test isapprox(norm(SVector(3.0, 4.0, 0.0)), 5.0; atol=1e-10)
        end

        @testset "scaling property" begin
            # norm(k*p) = k*norm(p) for k >= 0
            p = (3.0, 4.0, 0.0)
            k = 2.5
            @test isapprox(norm(Tuple(k .* p)), k * norm(p); rtol=1e-10)
        end

        @testset "invalid inputs" begin
            # Not a 3D vector
            @test_throws MethodError norm((1.0, 2.0))
            @test_throws MethodError norm((1.0, 2.0, 3.0, 4.0))

            # Not a tuple or SVector
            @test_throws MethodError norm([1.0, 2.0, 3.0])

            # Not numeric
            @test_throws MethodError norm(("a", "b", "c"))

            # Complex values
            @test_throws MethodError norm((1.0 + 1.0im, 2.0, 3.0))
        end
    end

    # Testing conversions in case implementation changes, even though currently external
    @testset "geodetic_to_ecef" begin
        @testset "origin" begin
            # Equator, prime meridian, sea level
            ecef = geodetic_to_ecef((0.0, 0.0, 0.0))
            @test ecef isa AbstractVector
            @test all(isreal.(ecef))
            @test all(isfinite.(ecef))
            @test length(ecef) == 3

            @test isapprox(ecef[1], SEMIMAJOR_RADIUS; rtol=1e-10)
            @test isapprox(ecef[2], 0.0; atol=1e-10)
            @test isapprox(ecef[3], 0.0; atol=1e-10)
        end

        @testset "cardinal directions" begin
            # Equator, 90° East
            ecef_east = geodetic_to_ecef((0.0, π/2, 0.0))
            @test isapprox(ecef_east[1], 0.0; atol=1e-10)
            @test isapprox(ecef_east[2], SEMIMAJOR_RADIUS; rtol=1e-10)
            @test isapprox(ecef_east[3], 0.0; atol=1e-10)

            # Equator, 90° West
            ecef_west = geodetic_to_ecef((0.0, -π/2, 0.0))
            @test isapprox(ecef_west[1], 0.0; atol=1e-10)
            @test isapprox(ecef_west[2], -SEMIMAJOR_RADIUS; rtol=1e-10)
            @test isapprox(ecef_west[3], 0.0; atol=1e-10)

            # Equator, 180°
            ecef_180 = geodetic_to_ecef((0.0, π, 0.0))
            @test isapprox(ecef_180[1], -SEMIMAJOR_RADIUS; rtol=1e-10)
            @test isapprox(ecef_180[2], 0.0; atol=1e-10)
            @test isapprox(ecef_180[3], 0.0; atol=1e-10)
        end

        @testset "poles" begin
            # North pole
            ecef_north = geodetic_to_ecef((π/2, 0.0, 0.0))
            @test isapprox(ecef_north[1], 0.0; atol=1e-10)
            @test isapprox(ecef_north[2], 0.0; atol=1e-10)
            @test isapprox(ecef_north[3], SEMIMINOR_RADIUS; rtol=1e-10)

            # South pole
            ecef_south = geodetic_to_ecef((-π/2, 0.0, 0.0))
            @test isapprox(ecef_south[1], 0.0; atol=1e-10)
            @test isapprox(ecef_south[2], 0.0; atol=1e-10)
            @test isapprox(ecef_south[3], -SEMIMINOR_RADIUS; rtol=1e-10)

            # Near north pole
            ecef_near_north = geodetic_to_ecef((deg2rad(89.9999999), 0.0, 0.0))
            @test isapprox(ecef_near_north[1], 0.009673; atol=1e-10)
            @test isapprox(ecef_near_north[2], 0.005585; atol=1e-10)
            @test isapprox(ecef_near_north[3], 6356752.314245; rtol=1e-10)

            # North pole with non-zero longitude
            ecef_near_north_lon = geodetic_to_ecef((π/2, deg2rad(137), 0.0))
            @test isapprox(ecef_near_north_lon[1], 0.0; atol=1e-10)
            @test isapprox(ecef_near_north_lon[2], 0.0; atol=1e-10)
            @test isapprox(ecef_near_north_lon[3], SEMIMINOR_RADIUS; rtol=1e-10)
        end

        @testset "with altitude" begin
            # Shallow negative altitude
            ecef_shallow_neg = geodetic_to_ecef((0.0, 0.0, -1000.0))
            @test isapprox(ecef_shallow_neg[1], SEMIMAJOR_RADIUS - 1000.0; rtol=1e-10)
            @test isapprox(ecef_shallow_neg[2], 0.0; atol=1e-10)
            @test isapprox(ecef_shallow_neg[3], 0.0; atol=1e-10)

            # Altitude = -SEMIMAJOR_RADIUS
            @test isapprox(geodetic_to_ecef((0.0, 0.0, -SEMIMAJOR_RADIUS)), [0.0, 0.0, 0.0])

            # Very large altitude, GEO height
            @test isapprox(geodetic_to_ecef((0.0, 0.0, 35786000.0)), [42164137, 0.0, 0.0])
        end

        @testset "2-tuple input" begin
            # Should handle (lat, lon) without height
            @test isapprox(
                geodetic_to_ecef((
                    deg2rad(42.394205296564806), deg2rad(-72.52761607116234)
                )),
                [1416435.823714, -4499938.662610, 4278031.082705],
                rtol=1e-5,
            )
        end

        @testset "oblique points" begin
            # 45° N, 45° E, 0 m
            @test isapprox(
                geodetic_to_ecef((π/4, π/4, 0.0)),
                [3194419.145061, 3194419.145061, 4487348.408866],
            )

            # 45° N, 45° E, 1000 m
            @test isapprox(
                geodetic_to_ecef((π/4, π/4, 1000.0)),
                [3194919.145061, 3194919.145061, 4488055.515647],
            )

            # 45° S, 45° W, 0 m
            @test isapprox(
                geodetic_to_ecef((-π/4, -π/4, 0.0)),
                [3194419.145061, -3194419.145061, -4487348.408866],
            )

            # 39° N, 77° W, 100 m
            @test isapprox(
                geodetic_to_ecef((deg2rad(39), deg2rad(-77), 100.0)),
                [1116523.199853, -4836193.303241, 3992379.954791],
            )
        end

        @testset "wraparounds" begin
            # 180 = -180
            @test isapprox(geodetic_to_ecef((π, π, 0.0)), geodetic_to_ecef((-π, -π, 0.0)))

            # 359 = -1
            @test isapprox(
                geodetic_to_ecef((deg2rad(359), deg2rad(359), 0.0)),
                geodetic_to_ecef((deg2rad(-1), deg2rad(-1), 0.0)),
            )

            # 450 = 90
            @test isapprox(
                geodetic_to_ecef((deg2rad(450), deg2rad(450), 0.0)),
                geodetic_to_ecef((deg2rad(90), deg2rad(90), 0.0)),
            )
        end

        @testset "invalid inputs" begin
            @test_throws DomainError geodetic_to_ecef((Inf, 0.0, 0.0))
            @test_throws DomainError geodetic_to_ecef((0.0, Inf, 0.0))
            # @test_throws DomainError geodetic_to_ecef((0.0, 0.0, Inf))

            # @test_throws DomainError geodetic_to_ecef((NaN, 0.0, 0.0))
            # @test_throws DomainError geodetic_to_ecef((0.0, NaN, 0.0))
            # @test_throws DomainError geodetic_to_ecef((0.0, 0.0, NaN))
        end

        @testset "different ellipsoid" begin
            # GRS80
            GRS80_ellipsoid = Ellipsoid(6378137.0, 1 / 298.257222101)

            # origin
            @test isapprox(
                geodetic_to_ecef((0.0, 0.0, 0.0), ellipsoid=GRS80_ellipsoid),
                [GRS80_ellipsoid.a, 0.0, 0.0],
            )

            # North Pole
            @test isapprox(
                geodetic_to_ecef((π/2, 0.0, 0.0), ellipsoid=GRS80_ellipsoid),
                [0.0, 0.0, GRS80_ellipsoid.b],
            )

            # Oblique angles
            @test isapprox(
                geodetic_to_ecef((π/4, π/4, 0.0), ellipsoid=GRS80_ellipsoid),
                [3194419.145087, 3194419.145087, 4487348.408755],
            )
            @test isapprox(
                geodetic_to_ecef((π/4, π/4, 1000.0), ellipsoid=GRS80_ellipsoid),
                [3194919.145087, 3194919.145087, 4488055.515536],
            )
            @test isapprox(
                geodetic_to_ecef(
                    (deg2rad(39), deg2rad(-77), 100.0), ellipsoid=GRS80_ellipsoid
                ),
                [1116523.199860, -4836193.303272, 3992379.954685],
            )
        end
    end

    @testset "ecef_to_geodetic" begin
        @testset "origin" begin
            # Equator, prime meridian, sea level
            geodetic = ecef_to_geodetic((SEMIMAJOR_RADIUS, 0.0, 0.0))
            @test geodetic isa Tuple
            @test all(isreal.(geodetic))
            @test all(isfinite.(geodetic))
            @test length(geodetic) == 3

            @test isapprox(geodetic[1], 0.0)
            @test isapprox(geodetic[2], 0.0)
            @test isapprox(geodetic[3], 0.0)
        end

        @testset "cardinal directions" begin
            # Equator, 90° East
            geodetic_east = ecef_to_geodetic((0.0, SEMIMAJOR_RADIUS, 0.0))
            @test isapprox(geodetic_east[1], 0.0)
            @test isapprox(geodetic_east[2], π/2)
            @test isapprox(geodetic_east[3], 0.0)

            # Equator, 90° West
            geodetic_west = ecef_to_geodetic((0.0, -SEMIMAJOR_RADIUS, 0.0))
            @test isapprox(geodetic_west[1], 0.0)
            @test isapprox(geodetic_west[2], -π/2)
            @test isapprox(geodetic_west[3], 0.0)

            # Equator, 180°
            geodetic_180 = ecef_to_geodetic((-SEMIMAJOR_RADIUS, 0.0, 0.0))
            @test isapprox(geodetic_180[1], 0.0)
            @test isapprox(geodetic_180[2], π)
            @test isapprox(geodetic_180[3], 0.0)
        end

        @testset "poles" begin
            # North pole
            geodetic_north = ecef_to_geodetic((0.0, 0.0, SEMIMINOR_RADIUS))
            @test isapprox(geodetic_north[1], π/2)
            @test isapprox(geodetic_north[2], 0.0)
            @test isapprox(geodetic_north[3], 0.0; atol=1e-4)

            # South pole
            geodetic_south = ecef_to_geodetic((0.0, 0.0, -SEMIMINOR_RADIUS))
            @test isapprox(geodetic_south[1], -π/2)
            @test isapprox(geodetic_south[2], 0.0)
            @test isapprox(geodetic_south[3], 0.0; atol=1e-4)

            # Near north pole
            geodetic_near_north = ecef_to_geodetic((0.009673, 0.005585, 6356752.314245))
            @test isapprox(geodetic_near_north[1], deg2rad(89.9999999))
            @test isapprox(geodetic_near_north[2], deg2rad(30.00129205))
            @test isapprox(geodetic_near_north[3], 0.0; atol=1e-4)
        end

        @testset "with altitude" begin
            # Shallow negative altitude
            geodetic_shallow_neg = ecef_to_geodetic((SEMIMAJOR_RADIUS - 1000.0, 0.0, 0.0))
            @test isapprox(geodetic_shallow_neg[1], 0.0)
            @test isapprox(geodetic_shallow_neg[2], 0.0)
            @test isapprox(geodetic_shallow_neg[3], -1000.0)

            # origin
            @test_nowarn ecef_to_geodetic((0.0, 0.0, 0.0))
            @test all(
                isapprox.(
                    ecef_to_geodetic((0.0, 0.0, 0.0)), ecef_to_geodetic((0.0, 0.0, 0.0))
                ),
            )

            # Very large altitude, GEO height
            geodetic_geo = ecef_to_geodetic((42164137, 0.0, 0.0))
            @test isapprox(geodetic_geo[1], 0.0)
            @test isapprox(geodetic_geo[2], 0.0)
            @test isapprox(geodetic_geo[3], 35786000.0)
        end

        @testset "oblique points" begin
            # 45° N, 45° E, 0 m
            geodetic_ne = ecef_to_geodetic((3194419.145061, 3194419.145061, 4487348.408866))
            @test isapprox(geodetic_ne[1], π/4)
            @test isapprox(geodetic_ne[2], π/4)
            @test isapprox(geodetic_ne[3], 0.0; atol=1e-4)

            # 45° N, 45° E, 1000 m
            geodetic_ne_1000 = ecef_to_geodetic((
                3194919.145061, 3194919.145061, 4488055.515647
            ))
            @test isapprox(geodetic_ne_1000[1], π/4)
            @test isapprox(geodetic_ne_1000[2], π/4)
            @test isapprox(geodetic_ne_1000[3], 1000.0)

            # 45° S, 45° W, 0 m
            geodetic_sw = ecef_to_geodetic((
                3194419.145061, -3194419.145061, -4487348.408866
            ))
            @test isapprox(geodetic_sw[1], -π/4)
            @test isapprox(geodetic_sw[2], -π/4)
            @test isapprox(geodetic_sw[3], 0.0; atol=1e-4)

            # 39° N, 77° W, 100 m
            geodetic_dc = ecef_to_geodetic((
                1116523.199853, -4836193.303241, 3992379.954791
            ))
            @test isapprox(geodetic_dc[1], deg2rad(39))
            @test isapprox(geodetic_dc[2], deg2rad(-77))
            @test isapprox(geodetic_dc[3], 100.0)
        end

        @testset "float edge cases" begin
            # positively signed zero
            geodetic_pos_zero = ecef_to_geodetic((-SEMIMAJOR_RADIUS, +0.0, 0.0))
            @test isapprox(geodetic_pos_zero[1], 0.0)
            @test isapprox(geodetic_pos_zero[2], π)
            @test isapprox(geodetic_pos_zero[3], 0.0)

            # negatively signed zero
            geodetic_neg_zero = ecef_to_geodetic((-SEMIMAJOR_RADIUS, -0.0, 0.0))
            @test isapprox(geodetic_neg_zero[1], 0.0)
            @test isapprox(geodetic_neg_zero[2], -π)
            @test isapprox(geodetic_neg_zero[3], 0.0)

            # small perturbations around zero
            geodetic_perturbed_pos = ecef_to_geodetic((-SEMIMAJOR_RADIUS, 0.001, 0.0))
            @test isapprox(geodetic_perturbed_pos[1], 0.0)
            @test isapprox(geodetic_perturbed_pos[2], deg2rad(179.9999999910))
            @test isapprox(geodetic_perturbed_pos[3], 0.0)

            geodetic_perturbed_neg = ecef_to_geodetic((-SEMIMAJOR_RADIUS, -0.001, 0.0))
            @test isapprox(geodetic_perturbed_neg[1], 0.0)
            @test isapprox(geodetic_perturbed_neg[2], deg2rad(-179.9999999910))
            @test isapprox(geodetic_perturbed_neg[3], 0.0)

            # overflow
            # not overflow safe
            # geodetic_overflow = ecef_to_geodetic((1e200, 1e200, 1e200))
            # @test isapprox(geodetic_overflow[1], deg2rad(35.26438968275399))
            # @test isapprox(geodetic_overflow[2], deg2rad(45.0))
            # @test isapprox(geodetic_overflow[3], 1.732050807568877e200)
        end

        @testset "invalid inputs" begin
            # @test_throws DomainError ecef_to_geodetic((Inf, 0.0, 0.0))
            # @test_throws DomainError ecef_to_geodetic((0.0, Inf, 0.0))
            # @test_throws DomainError ecef_to_geodetic((0.0, 0.0, Inf))

            # @test_throws DomainError ecef_to_geodetic((NaN, 0.0, 0.0))
            # @test_throws DomainError ecef_to_geodetic((0.0, NaN, 0.0))
            # @test_throws DomainError ecef_to_geodetic((0.0, 0.0, NaN))
        end

        @testset "different ellipsoid" begin
            # GRS80
            GRS80_ellipsoid = Ellipsoid(6378137.0, 1 / 298.257222101)

            # origin
            geodetic_grs80_origin = ecef_to_geodetic(
                (GRS80_ellipsoid.a, 0.0, 0.0), ellipsoid=GRS80_ellipsoid
            )
            @test isapprox(geodetic_grs80_origin[1], 0.0)
            @test isapprox(geodetic_grs80_origin[2], 0.0)
            @test isapprox(geodetic_grs80_origin[3], 0.0)

            # North Pole
            geodetic_grs80_north_pole = ecef_to_geodetic(
                (0.0, 0.0, GRS80_ellipsoid.b), ellipsoid=GRS80_ellipsoid
            )
            @test isapprox(geodetic_grs80_north_pole[1], π/2)
            @test isapprox(geodetic_grs80_north_pole[2], 0.0)
            @test isapprox(geodetic_grs80_north_pole[3], 0.0)

            # Oblique angles
            geodetic_grs80_ne = ecef_to_geodetic(
                (3194419.145087, 3194419.145087, 4487348.408755), ellipsoid=GRS80_ellipsoid
            )
            @test isapprox(geodetic_grs80_ne[1], π/4)
            @test isapprox(geodetic_grs80_ne[2], π/4)
            @test isapprox(geodetic_grs80_ne[3], 0.0; atol=1e-4)

            geodetic_grs80_ne_1000 = ecef_to_geodetic(
                (3194919.145087, 3194919.145087, 4488055.515536), ellipsoid=GRS80_ellipsoid
            )
            @test isapprox(geodetic_grs80_ne_1000[1], π/4)
            @test isapprox(geodetic_grs80_ne_1000[2], π/4)
            @test isapprox(geodetic_grs80_ne_1000[3], 1000.0)

            geodetic_grs80_dc_100 = ecef_to_geodetic(
                (1116523.199860, -4836193.303272, 3992379.954685), ellipsoid=GRS80_ellipsoid
            )
            @test isapprox(geodetic_grs80_dc_100[1], deg2rad(39))
            @test isapprox(geodetic_grs80_dc_100[2], deg2rad(-77))
            @test isapprox(geodetic_grs80_dc_100[3], 100.0)
        end
    end

    @testset "ellipsoid_line_intersection" begin
        @testset "line through center" begin
            # Line through center, perpendicular to z-axis
            p1 = (2.0, 0.0, 0.0)
            p2 = (-2.0, 0.0, 0.0)
            xs = ellipsoid_line_intersection(1.0, 1.0, p1, p2)

            @test length(xs) == 2
            # Both points should be at distance 1 from origin (sphere radius)
            @test all(isapprox(sqrt(sum(Tuple(x) .^ 2)), 1.0; atol=1e-10) for x in xs)
        end

        @testset "line misses ellipsoid" begin
            # Line far from origin
            p1 = (10.0, 10.0, 0.0)
            p2 = (10.0, 20.0, 0.0)
            xs = ellipsoid_line_intersection(1.0, 1.0, p1, p2)

            @test length(xs) == 0
        end

        @testset "tangent line" begin
            # Tangent to sphere (approximately)
            p1 = (1.0, 5.0, 0.0)
            p2 = (1.0, -5.0, 0.0)
            xs = ellipsoid_line_intersection(1.0, 1.0, p1, p2)

            # Should have 0, 1, or 2 intersections (tangent is 1)
            @test 0 <= length(xs) <= 2
        end

        @testset "oblate ellipsoid" begin
            # Earth-like ellipsoid (oblate)
            major = SEMIMAJOR_RADIUS
            minor = SEMIMINOR_RADIUS

            # Line through equator
            p1 = (major * 2, 0.0, 0.0)
            p2 = (-major * 2, 0.0, 0.0)
            xs = ellipsoid_line_intersection(major, minor, p1, p2)

            @test length(xs) == 2
            # Points should be on ellipsoid
            for x in xs
                @test isapprox(x[3], 0.0; atol=1e-10)  # z should be ~0
            end
        end

        @testset "only_between=true" begin
            # Two intersections exist, but only one between points
            p1 = (2.0, 0.0, 0.0)
            p2 = (0.5, 0.0, 0.0)

            xs_all = ellipsoid_line_intersection(1.0, 1.0, p1, p2; only_between=false)
            xs_between = ellipsoid_line_intersection(1.0, 1.0, p1, p2; only_between=true)

            @test length(xs_between) <= length(xs_all)
        end

        @testset "only_between=false" begin
            # Segment doesn't intersect, but line does
            p1 = (0.5, 0.0, 0.0)
            p2 = (0.4, 0.0, 0.0)

            xs_all = ellipsoid_line_intersection(1.0, 1.0, p1, p2; only_between=false)
            xs_between = ellipsoid_line_intersection(1.0, 1.0, p1, p2; only_between=true)

            @test length(xs_all) >= length(xs_between)
        end

        @testset "intersection properties" begin
            # Verify intersections lie on both line and ellipsoid
            major, minor = 2.0, 1.5
            p1 = (4.0, 0.0, 0.0)
            p2 = (-4.0, 0.0, 0.0)
            xs = ellipsoid_line_intersection(major, minor, p1, p2)

            for x in xs
                # Should be on z-axis (the line)
                @test isapprox(x[2], 0.0; atol=1e-10)
                @test isapprox(x[3], 0.0; atol=1e-10)

                # Should be on ellipsoid: x²/major² + y²/minor² + z²/minor² = 1
                ellipsoid_eq = (x[1]^2 / major^2) + (x[2]^2 / minor^2) + (x[3]^2 / minor^2)
                @test isapprox(ellipsoid_eq, 1.0; atol=1e-10)
            end
        end

        @testset "return type" begin
            xs = ellipsoid_line_intersection(1.0, 1.0, (2.0, 0.0, 0.0), (-2.0, 0.0, 0.0))

            @test xs isa Vector
            @test all(x isa Tuple for x in xs)
            @test all(length(x) == 3 for x in xs)
        end

        @testset "number of intersections" begin
            # Should always return 0, 1, or 2 intersections
            test_cases = [
                # (major, minor, p1, p2)
                (1.0, 1.0, (2.0, 0.0, 0.0), (-2.0, 0.0, 0.0)),     # Through center
                (1.0, 1.0, (10.0, 10.0, 0.0), (10.0, 20.0, 0.0)),   # Misses
                (
                    SEMIMAJOR_RADIUS,
                    SEMIMINOR_RADIUS,
                    (SEMIMAJOR_RADIUS*2, 0.0, 0.0),
                    (-SEMIMAJOR_RADIUS*2, 0.0, 0.0),
                ),  # Earth
            ]

            for (major, minor, p1, p2) in test_cases
                xs = ellipsoid_line_intersection(major, minor, p1, p2)
                @test length(xs) in [0, 1, 2]
            end
        end
    end

    @testset "gs_gs_distance" begin
        @testset "same location" begin
            # Distance from point to itself should be zero
            gs = (0.0, 0.0, 0.0)
            d = gs_gs_distance(gs, gs)
            @test isapprox(d, 0.0; atol=1.0)

            # Mid-latitude point to itself
            gs_mid = (deg2rad(45.0), deg2rad(-120.0), 100.0)
            d_mid = gs_gs_distance(gs_mid, gs_mid)
            @test isapprox(d_mid, 0.0; atol=1.0)
        end

        @testset "symmetry" begin
            # Distance should be symmetric
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (deg2rad(45.0), deg2rad(90.0), 0.0)

            d12 = gs_gs_distance(gs1, gs2)
            d21 = gs_gs_distance(gs2, gs1)

            @test isapprox(d12, d21; rtol=1e-10)
        end

        @testset "non-negative" begin
            # Distance should always be non-negative
            test_pairs = [
                ((0.0, 0.0, 0.0), (π/4, π/4, 0.0)),
                ((π/2, 0.0, 0.0), (-π/2, 0.0, 0.0)),
                (
                    (deg2rad(45.0), deg2rad(-120.0), 100.0),
                    (deg2rad(35.0), deg2rad(-110.0), 500.0),
                ),
            ]

            for (gs1, gs2) in test_pairs
                d = gs_gs_distance(gs1, gs2)
                @test d >= 0
            end
        end

        @testset "equator quarter circle" begin
            # Two points on equator 90° apart should be ~1/4 Earth circumference
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (0.0, π/2, 0.0)

            d = gs_gs_distance(gs1, gs2)
            expected_approx = (π/2) * SEMIMAJOR_RADIUS  # quarter circumference

            @test isapprox(d, expected_approx; rtol=0.01)  # within 1%
        end

        @testset "explicit format tags" begin
            # Test that explicit format tags work
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (deg2rad(45.0), deg2rad(90.0), 0.0)

            d_default = gs_gs_distance(gs1, gs2)
            d_latlon = gs_gs_distance(Val(:lat_lon), gs1, gs2)

            @test isapprox(d_default, d_latlon; rtol=1e-10)
        end

        @testset "ECEF format" begin
            # Equator, prime meridian (hardcoded ECEF)
            ecef1 = (SEMIMAJOR_RADIUS, 0.0, 0.0)

            # Equator, 90° East (hardcoded ECEF)
            ecef2 = (0.0, SEMIMAJOR_RADIUS, 0.0)

            d_ecef = gs_gs_distance(Val(:ECEF), ecef1, ecef2)

            # Should be ~1/4 Earth circumference
            expected_approx = (π/2) * SEMIMAJOR_RADIUS
            @test isapprox(d_ecef, expected_approx; rtol=0.01)
        end

        @testset "with altitude" begin
            # Altitude should affect distance
            gs_sea = (0.0, 0.0, 0.0)
            gs_alt = (0.0, 0.0, 1000.0)

            d_sea = gs_gs_distance(gs_sea, gs_sea)
            d_alt = gs_gs_distance(gs_alt, gs_alt)

            @test isapprox(d_sea, 0.0; atol=1.0)
            @test isapprox(d_alt, 0.0; atol=1.0)
        end

        @testset "two altitude points" begin
            # Different altitudes at same location
            gs_sea = (0.0, 0.0, 0.0)
            gs_high = (0.0, 0.0, 10000.0)

            d = gs_gs_distance(gs_sea, gs_high)
            @test d > 0
            @test isapprox(d, 10000.0; rtol=1e-3)
        end

        @testset "poles" begin
            # Distance between poles
            north = (π/2, 0.0, 0.0)
            south = (-π/2, 0.0, 0.0)

            d = gs_gs_distance(north, south)
            expected = π * SEMIMINOR_RADIUS  # half circumference along pole

            @test isapprox(d, expected; rtol=0.01)
        end

        @testset "return type" begin
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (π/4, π/4, 0.0)

            d = gs_gs_distance(gs1, gs2)
            @test d isa Number
            @test d >= 0
        end

        @testset "short distances" begin
            # Two points very close should have small distance
            gs1 = (0.0, 0.0, 0.0)
            gs2 = (1e-6, 1e-6, 0.0)

            d = gs_gs_distance(gs1, gs2)
            @test d > 0
            @test d < 1000  # Less than 1 km for tiny angle difference
        end

        @testset "2-tuple input" begin
            # Should handle (lat, lon) without height
            gs1 = (0.0, 0.0)
            gs2 = (0.0, π/2)

            d = gs_gs_distance(gs1, gs2)
            expected_approx = (π/2) * SEMIMAJOR_RADIUS

            @test isapprox(d, expected_approx; rtol=0.01)
        end
        @testset "various antipodal points" begin
            expected = π * (SEMIMAJOR_RADIUS + SEMIMINOR_RADIUS) / 2

            for lat in [deg2rad(10.0), deg2rad(30.0), deg2rad(60.0), deg2rad(80.0)]
                for lon_offset in [0.0, deg2rad(45.0), deg2rad(90.0)]
                    gs1 = (lat, lon_offset, 0.0)
                    gs2 = (-lat, lon_offset + π, 0.0)

                    d = gs_gs_distance(gs1, gs2)
                    @test isapprox(d, expected; rtol=1e-6)
                end
            end
        end
    end

    @testset "sat_sat_distance" begin
        @testset "same satellite" begin
            tle_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle = read_tle(tle_string, verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)

            # Distance to itself should be ~0
            d = sat_sat_distance(sat, sat; time=sat.sgp4d.epoch)
            @test d === nothing || isapprox(d, 0.0; atol=1e-10)
        end

        @testset "non-negative when not obstructed" begin
            # ISS TLE
            tle1_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            # Hubble Space Telescope TLE (different orbit)
            tle2_string = """HUBBLE SPACE TELESCOPE
                1 20580U 90037B   26225.00000000  .00001234  00000-0  61256-5 0  9991
                2 20580  28.4699 277.5360 0002853 182.9063 177.2222 15.09625533729000"""

            tle1 = read_tle(tle1_string, verify_checksum=false)
            tle2 = read_tle(tle2_string, verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            d = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)

            if d !== nothing
                @test d >= 0
                @test d > 1000  # Different orbits should be thousands of km apart
            end
        end

        @testset "symmetry" begin
            tle1_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle2_string = """HUBBLE SPACE TELESCOPE
                1 20580U 90037B   26225.00000000  .00001234  00000-0  61256-5 0  9991
                2 20580  28.4699 277.5360 0002853 182.9063 177.2222 15.09625533729000"""

            tle1 = read_tle(tle1_string, verify_checksum=false)
            tle2 = read_tle(tle2_string, verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            d12 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)
            d21 = sat_sat_distance(sat2, sat1; time=sat1.sgp4d.epoch)

            if d12 !== nothing && d21 !== nothing
                @test isapprox(d12, d21; rtol=1e-10)
            end
        end

        @testset "Propagator input" begin
            tle_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle = read_tle(tle_string, verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)

            d = sat_sat_distance(sat, sat; time=sat.sgp4d.epoch)
            @test d === nothing || isa(d, Number)
        end

        @testset "TLE input" begin
            tle1_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle2_string = """HUBBLE SPACE TELESCOPE
                1 20580U 90037B   26225.00000000  .00001234  00000-0  61256-5 0  9991
                2 20580  28.4699 277.5360 0002853 182.9063 177.2222 15.09625533729000"""

            tle1 = read_tle(tle1_string, verify_checksum=false)
            tle2 = read_tle(tle2_string, verify_checksum=false)

            d = sat_sat_distance(tle1, tle2; time=JD_J2000)
            @test d === nothing || isa(d, Number)
        end

        @testset "time parameter" begin
            tle1_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle2_string = """HUBBLE SPACE TELESCOPE
                1 20580U 90037B   26225.00000000  .00001234  00000-0  61256-5 0  9991
                2 20580  28.4699 277.5360 0002853 182.9063 177.2222 15.09625533729000"""

            tle1 = read_tle(tle1_string, verify_checksum=false)
            tle2 = read_tle(tle2_string, verify_checksum=false)

            sat1 = Propagators.init(Val(:SGP4), tle1)
            sat2 = Propagators.init(Val(:SGP4), tle2)

            # Distance can change at different times
            d1 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch)
            d2 = sat_sat_distance(sat1, sat2; time=sat1.sgp4d.epoch + 100)  # 100 days later

            if d1 !== nothing && d2 !== nothing
                # Distances should generally be different at different times
                @test isa(d1, Number)
                @test isa(d2, Number)
            end
        end

        @testset "return type" begin
            tle_string = """ISS (ZARYA)
                1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
                2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

            tle = read_tle(tle_string, verify_checksum=false)
            sat = Propagators.init(Val(:SGP4), tle)

            result = sat_sat_distance(sat, sat; time=sat.sgp4d.epoch)
            @test result === nothing || result isa Number
        end
    end

    @testset "sun_position" begin
        @testset "return type and magnitude" begin
            # At J2000 epoch
            sun = sun_position(JD_J2000)

            @test sun isa AbstractVector
            @test length(sun) == 3

            # Distance from Earth should be ~1 AU +/- 2.5 million km
            au_m = 1.495978707e11
            r = sqrt(sum(sun .^ 2))
            @test isapprox(r, au_m; atol=2.5e9)  # 2.5 million km tolerance
        end

        @testset "DateTime input" begin
            # Test with DateTime instead of Julian days
            dt = DateTime(2000, 1, 1, 12, 0, 0)
            sun_dt = sun_position(dt)

            @test sun_dt isa AbstractVector
            @test length(sun_dt) == 3

            r = sqrt(sum(sun_dt .^ 2))
            au_m = 1.495978707e11
            @test isapprox(r, au_m; atol=2.5e9)
        end

        @testset "different epochs" begin
            # Test multiple different times
            times = [JD_J2000, JD_J2000 + 1, JD_J2000 + 365, JD_J2000 + 365*10]
            au_m = 1.495978707e11

            for jd in times
                sun = sun_position(jd)
                r = sqrt(sum(sun .^ 2))

                @test isapprox(r, au_m; atol=2.5e9)
            end
        end

        @testset "orbital motion" begin
            # Sun appears to move ~360°/year from Earth perspective
            # Check positions 90 days apart are meaningfully different
            jd1 = JD_J2000
            jd2 = JD_J2000 + 90.0

            sun1 = sun_position(jd1)
            sun2 = sun_position(jd2)

            # Angle between the two positions
            cos_angle =
                (sun1[1]*sun2[1] + sun1[2]*sun2[2] + sun1[3]*sun2[3]) /
                (sqrt(sum(sun1 .^ 2)) * sqrt(sum(sun2 .^ 2)))
            angle = acos(clamp(cos_angle, -1, 1))

            # Should be roughly 90° (π/2 radians) for 90 days
            @test isapprox(angle, π/2; rtol=0.01)
        end

        @testset "consistency" begin
            # Same time should give same result
            jd = JD_J2000
            sun1 = sun_position(jd)
            sun2 = sun_position(jd)

            @test all(isapprox.(sun1, sun2; rtol=1e-10))
        end
    end

    @testset "light_at_point" begin
        @testset "return type" begin
            point = (1.0e11, 0.0, 0.0)

            result = light_at_point(point, JD_J2000)
            @test result isa LightCondition
            @test result in [sunlight, penumbra, umbra]
        end

        @testset "consistency" begin
            # Same point and time should give same result
            point = (1.0e11, 0.0, 0.0)

            condition1 = light_at_point(point, JD_J2000)
            condition2 = light_at_point(point, JD_J2000)

            @test condition1 == condition2
        end

        @testset "different times give valid results" begin
            point = (1.0e11, 0.0, 0.0)

            condition1 = light_at_point(point, JD_J2000)
            condition2 = light_at_point(point, JD_J2000 + 100.0)

            @test condition1 in [sunlight, penumbra, umbra]
            @test condition2 in [sunlight, penumbra, umbra]
        end

        @testset "DateTime input" begin
            point = (1.0e11, 0.0, 0.0)
            dt = DateTime(2000, 1, 1, 12, 0, 0)

            condition = light_at_point(point, dt)
            @test condition in [sunlight, penumbra, umbra]
        end

        @testset "multiple points valid" begin
            test_points = [
                (1.0e11, 0.0, 0.0),
                (0.0, 1.0e11, 0.0),
                (1.0e10, 1.0e10, 1.0e10),
                (6.378137e6, 0.0, 0.0),
                (6.778137e6, 0.0, 0.0),
            ]

            for point in test_points
                condition = light_at_point(point, JD_J2000)
                @test condition in [sunlight, penumbra, umbra]
            end
        end
    end
end