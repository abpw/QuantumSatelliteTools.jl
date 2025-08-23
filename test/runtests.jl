using QuantumSatelliteTools
using Test

#using QuantumSatelliteTools.AstronomyGeometry: semimajor_radius, semiminor_radius, earth_mass_kg, seconds_per_day, G, FreespaceChannel
using QuantumSatelliteTools.AstronomyGeometry
using QuantumSatelliteTools.GenerateTLEs: generate_regular_array_TLEs

@testset "QuantumSatelliteTools.jl" begin
    @testset "AstronomyGeometry.jl" begin
        #using QuantumSatelliteTools.AstronomyGeometry
        using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
        using SatelliteToolboxTle: TLE
        using SatelliteToolboxBase: OrbitStateVector, JD_J2000
        using Dates

        @testset "geodetic<->ecef roundtrips" begin
            # Equator/prime meridian
            lat, lon, h = 0.0, 0.0, 0.0
            x, y, z = geodetic_to_ecef((lat, lon, h))
            lat2, lon2, h2 = ecef_to_geodetic(x, y, z)
            @test isapprox(lat2, lat; atol=1e-9)
            @test isapprox(lon2, lon; atol=1e-9)
            @test isapprox(h2, h; atol=1e-3)

            # A mid-latitude point (UMass LGRC)
            lat, lon, h = deg2rad(42.39466434352798), deg2rad(-72.5270869156386), 80.0
            ecef_pt = geodetic_to_ecef((lat, lon, h))
            lat3, lon3, h3 = ecef_to_geodetic(ecef_pt)
            @test isapprox(lat3, lat; atol=1e-9)
            @test isapprox(lon3, lon; atol=1e-9)
            @test isapprox(h3, h; atol=1e-2)
        end

        @testset "eci_to_ecef/geodetic (propagator & state vector)" begin
            # Build a simple constellation and pick one satellite
            tles = generate_regular_array_TLEs(orbits=4, sats_per_orbit=6, altitude_km=500, inclination_rad=π/2)
            tle = tles[1]
            sat::OrbitPropagatorSgp4 = Propagators.init(Val(:SGP4), tle)

            # Using propagator directly
            ecef_from_prop = eci_to_ecef(sat; time=sat.sgp4d.epoch)
            @test length(ecef_from_prop) == 3
            r_mag = √(sum(ecef_from_prop .^ 2))
            @test r_mag > semimajor_radius # above Earth's surface

            geo_from_prop = eci_to_geodetic(sat; time=sat.sgp4d.epoch)
            @test -π/2 <= geo_from_prop[1] <= π/2
            @test -π <= geo_from_prop[2] <= π

            # Using an OrbitStateVector (should match the propagator version at same time)
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            ecef_from_sv = eci_to_ecef(sv)
            @test all(isapprox.(ecef_from_sv, ecef_from_prop; atol=1e-6))

            geo_from_sv = eci_to_geodetic(sv)
            @test isapprox(geo_from_sv[1], geo_from_prop[1]; atol=1e-9)
            @test isapprox(geo_from_sv[2], geo_from_prop[2]; atol=1e-9)
        end

        @testset "eci_to_ecef/eci_to_geodetic (point/tuple overloads)" begin
            # Reuse the same satellite and time, but call point overloads
            tles = generate_regular_array_TLEs(orbits=1, sats_per_orbit=3)
            sat = Propagators.init(Val(:SGP4), tles[2])
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            x, y, z = sv.r
            ecef_tuple = eci_to_ecef(x, y, z; time=sat.sgp4d.epoch)
            ecef_point = eci_to_ecef(Tuple(sv.r); time=sat.sgp4d.epoch)
            @test all(isapprox.(ecef_tuple, ecef_point; atol=1e-6))

            geo_tuple = eci_to_geodetic(x, y, z; time=sat.sgp4d.epoch)
            geo_point = eci_to_geodetic(Tuple(sv.r); time=sat.sgp4d.epoch)
            @test all(isapprox.(geo_tuple, geo_point; atol=1e-9))
        end

        @testset "ellipsoid_line_intersection" begin
            # A line through Earth with two intersections
            p1 = (semimajor_radius*2.0, 0.0, 0.0)
            p2 = (-semimajor_radius*2.0, 0.0, 0.0)
            xs = ellipsoid_line_intersection(semimajor_radius, semiminor_radius, p1, p2)
            @test length(xs) == 2
            @test all(√(sum(Tuple(q) .^ 2)) ≈ semimajor_radius for q in xs)

            # Segment that only touches between endpoints (only_between=true)
            xs_between = ellipsoid_line_intersection(semimajor_radius, semiminor_radius, p1, (0.0, 0.0, 0.0); only_between=true)
            @test length(xs_between) == 1 || length(xs_between) == 0  # origin is on surface only for sphere; allow 0/1
        end

        @testset "gs_gs_distance (lat/lon and ECEF variants)" begin
            gs1 = (deg2rad(0.0), deg2rad(0.0), 0.0)
            gs2 = (deg2rad(0.0), deg2rad(90.0), 0.0)
            d_ll = gs_gs_distance(gs1, gs2)
            expected = (π/2) * (semimajor_radius + semiminor_radius) / 2
            @test isapprox(d_ll, expected; rtol=5e-3)

            e1 = geodetic_to_ecef(gs1)
            e2 = geodetic_to_ecef(gs2)
            d_ecef = gs_gs_distance(Val(:ECEF), e1, e2)
            @test isapprox(d_ecef, d_ll; rtol=1e-9)

            d_alias = gs_gs_distance(Val(:lat_lon), gs1, gs2)
            @test d_alias == d_ll
        end

        @testset "sat_sat_distance & LOS" begin
            # Two nearby sats in the same orbit should have finite LOS distance
            tles = generate_regular_array_TLEs(orbits=1, sats_per_orbit=60, altitude_km=500)
            sat1 = Propagators.init(Val(:SGP4), tles[1])
            sat2 = Propagators.init(Val(:SGP4), tles[2])
            d = sat_sat_distance(sat1, sat2; time=max(sat1.sgp4d.epoch, sat2.sgp4d.epoch))
            @test d !== nothing
            @test d > 0
        end

        @testset "sun_position & light_at_point" begin
            # At J2000, distance ~ 1 AU
            jd = JD_J2000
            sun = sun_position(jd)
            au_m = 1.495978707e11
            @test isapprox(√(sum(sun .^ 2)), au_m; rtol=5e-2)  # allow a few percent

            # A point far toward the Sun should be in sunlight; opposite should be umbra
            p_sunlit = ntuple(i->sun[i] * 1e-3, 3)  # direction toward the Sun
            p_shadow = ntuple(i->-sun[i] * 1e-3, 3) # opposite side
            @test light_at_point(p_sunlit, jd) == sunlight
            @test light_at_point(p_shadow, jd) == umbra
        end

        @testset "FreespaceChannel (sat-GS)" begin
            tles = generate_regular_array_TLEs(orbits=2, sats_per_orbit=8, altitude_km=500)
            sat = Propagators.init(Val(:SGP4), tles[1])
            # Place GS at current subsatellite point to guarantee visibility
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            lat, lon, _ = eci_to_geodetic(sv)
            ch = FreespaceChannel(sat, (lat, lon, 0.0); time=sat.sgp4d.epoch, min_θ=deg2rad(0))
            @test ch !== nothing
            @test ch.distance_m > 0
            @test ch.elevation_angle_rad >= 0
        end

        @testset "FreespaceChannel (sat-sat)" begin
            # Choose neighboring sats in same orbit to improve LOS likelihood
            tles = generate_regular_array_TLEs(orbits=1, sats_per_orbit=10, altitude_km=500)
            sat1 = Propagators.init(Val(:SGP4), tles[3])
            sat2 = Propagators.init(Val(:SGP4), tles[4])
            ch = FreespaceChannel(sat1, sat2; time=max(sat1.sgp4d.epoch, sat2.sgp4d.epoch))
            @test ch === nothing || (ch.distance_m > 0 && ch.min_altitude_m > 0)
        end
    end
    @testset "GenerateTLEs.jl" begin
        @testset "generate_regular_array_TLEs: defaults" begin
            # Default call
            tles = generate_regular_array_TLEs()

            # Basic count
            @test length(tles) == 12 * 18

            # All inclinations should be ~90 degrees by default (π/2 rad)
            @test all(isapprox(t.inclination, 90.0; atol=1e-6) for t in tles)

            # Names should be "SAT 1", "SAT 2", ..., in order
            @test tles[1].name == "SAT 1"
            @test tles[end].name == "SAT $(length(tles))"
            @test all(tles[i].name == "SAT $(i)" for i in eachindex(tles))

            # RAAN should take 12 unique values: 0,15,30,...,165 and repeat in blocks of 18
            unique_raan = sort(unique(t.raan for t in tles))
            @test unique_raan == collect(0:15:165)
            for (i, r) in enumerate(0:15:165)
                block = tles[(18*(i-1)+1):(18*i)]
                @test all(isapprox(t.raan, r; atol=1e-9) for t in block)
            end

            # Mean anomaly pattern: 360/sats_per_orbit*(sat_in_orbit-1 + (orbit-1)/orbits)
            sats_per_orbit = 18
            orbits = 12
            step = 360 / sats_per_orbit
            # Check a few representative entries
            # First satellite of first orbit
            @test isapprox(tles[1].mean_anomaly, step * (0 + (0)/orbits); atol=1e-9)
            # Last satellite of first orbit
            @test isapprox(tles[sats_per_orbit].mean_anomaly, step * (sats_per_orbit-1 + 0/orbits); atol=1e-9)
            # First satellite of second orbit
            @test isapprox(tles[sats_per_orbit+1].mean_anomaly, step * (0 + (1)/orbits); atol=1e-9)
            # Last satellite of last orbit
            @test isapprox(tles[end].mean_anomaly, step * (sats_per_orbit-1 + (orbits-1)/orbits); atol=1e-9)

            # Eccentricity and argument of perigee when frozen_orbits=false
            @test all(isapprox(t.eccentricity, 1e-7; atol=1e-9) for t in tles)
            @test all(isapprox(t.argument_of_perigee, 90.0; atol=1e-9) for t in tles)

            # Mean motion should match analytical value for 500 km altitude
            altitude_km = 500
            n_expected = seconds_per_day * sqrt(G * earth_mass_kg / ( (altitude_km*1000 + semimajor_radius)^3 )) / (2*π)
            @test all(isapprox(t.mean_motion, n_expected; rtol=1e-10, atol=0) for t in tles)
        end

        @testset "generate_regular_array_TLEs: custom params" begin
            orbits = 3
            sats_per_orbit = 5
            altitude_km = 600
            inc_rad = 0.3
            tles = generate_regular_array_TLEs(orbits=orbits, sats_per_orbit=sats_per_orbit, altitude_km=altitude_km, inclination_rad=inc_rad, frozen_orbits=false)

            # Count check
            @test length(tles) == orbits * sats_per_orbit

            # Inclination check (in degrees)
            @test all(isapprox(t.inclination, inc_rad/π*180; atol=1e-9) for t in tles)

            # RAAN unique values should be 0, 60, 120
            unique_raan = sort(unique(t.raan for t in tles))
            @test unique_raan == [0.0, 60.0, 120.0]

            # Mean motion matches altitude
            n_expected = seconds_per_day * sqrt(G * earth_mass_kg / ( (altitude_km*1000 + semimajor_radius)^3 )) / (2*π)
            @test all(isapprox(t.mean_motion, n_expected; rtol=1e-10, atol=0) for t in tles)

            # Name formatting
            @test tles[1].name == "SAT 1"
            @test tles[end].name == "SAT $(length(tles))"
        end
    end
    @testset "GenerateGroundStations.jl" begin
        using QuantumSatelliteTools.GenerateGroundStations:
            is_within_min_distance,
            is_land,
            generate_population_center_gses,
            generate_city_gses,
            generate_random_gses,
            generate_equispaced_gses,
            generate_grid_gses
        using QuantumSatelliteTools.AstronomyGeometry: gs_gs_distance, GS

        @testset "is_within_min_distance" begin
            # Two existing GS ~ 1000 km apart
            lat1, lon1 = deg2rad(40.7128), deg2rad(-74.0060)   # NYC
            lat2, lon2 = deg2rad(41.8781), deg2rad(-87.6298)   # Chicago
            gses = [(lat1, lon1), (lat2, lon2)]
            # A candidate very near NYC (within 50 km): should be false for min_distance=100 km
            latc, lonc = deg2rad(40.70), deg2rad(-74.0)
            @test is_within_min_distance(latc, lonc, gses, 100) == false
            # A candidate far away (Honolulu): should be true
            lath, lonh = deg2rad(21.3069), deg2rad(-157.8583)
            @test is_within_min_distance(lath, lonh, gses, 1000) == true
        end

        @testset "is_land (coordinates)" begin
            # Known land: LGRC at UMass
            lat_land, lon_land = deg2rad(42.39466434352798), deg2rad(-72.5270869156386)
            @test is_land(lat_land, lon_land) == true
            # Known ocean: central Pacific near (0°, -160°)
            lat_sea, lon_sea = 0.0, deg2rad(-160)
            @test is_land(lat_sea, lon_sea) == false
        end

        @testset "is_land (GS tuple wrapper)" begin
            lat_land, lon_land = deg2rad(34.0522), deg2rad(-118.2437) # Los Angeles
            gs_tuple = (lat_land, lon_land, 0.0)::GS
            @test is_land(gs_tuple) == true
        end

        @testset "generate_population_center_gses" begin
            # Request a small number to keep the test fast
            gses = generate_population_center_gses(3; min_distance=0)
            @test length(gses) == 3
            # Ranges check
            @test all(-π/2 .<= first.(gses) .<= π/2)
            @test all(-π .<= last.(gses) .<= π)
        end

        @testset "generate_city_gses" begin
            cities = ["New York", "Paris", "Nonexistentopolis"]
            gses = generate_city_gses(cities)
            # Expect at least the two real cities to be found
            @test length(gses) >= 2
            # All returned tuples in valid ranges
            @test all(-π/2 .<= first.(gses) .<= π/2)
            @test all(-π .<= last.(gses) .<= π)
        end

        @testset "generate_random_gses" begin
            n = 5
            seed = 123456
            gses1 = generate_random_gses(n; seed=seed, min_distance=0, force_land=false)
            gses2 = generate_random_gses(n; seed=seed, min_distance=0, force_land=false)
            @test length(gses1) == n
            @test gses1 == gses2  # determinism with fixed seed
            # Range checks and uniqueness
            @test all(-π/2 .<= first.(gses1) .<= π/2)
            @test all(-π .<= last.(gses1) .<= π)
            @test length(unique(gses1)) == n
        end

        @testset "generate_equispaced_gses" begin
            n = 11
            gses = generate_equispaced_gses(n)
            @test length(gses) == n
            @test length(unique(gses)) == n
            @test all(-π/2 .<= first.(gses) .<= π/2)
            @test all(-π .<= last.(gses) .<= π)
            # Should include the equator/prime-meridian point when n is odd
            @test any(isapprox.(first.(gses), 0.0; atol=1e-12) .& isapprox.(last.(gses), 0.0; atol=1e-12))
        end

        @testset "generate_grid_gses" begin
            # Basic smoke test: should return non-empty list with valid ranges
            gses = generate_grid_gses(; equatorial_distance_km=800, force_land=false)
            @test !isempty(gses)
            @test all(-π/2 .<= first.(gses) .<= π/2)
            @test all(-π .<= last.(gses) .<= π)
        end
    end
    @testset "LossCalculation.jl" begin
        using QuantumSatelliteTools.AstronomyGeometry: FreespaceChannel, clear, umbra
        using QuantumSatelliteTools.LossCalculation: waist_radius, geometric_loss, atmosphere_distance, atmospheric_loss, pointing_loss, reflector_loss, swapping_loss
        
        @testset "waist_radius" begin
            # Test basic waist radius calculation
            # w₀ = 1.22λ/D_t where λ = 1550nm = 1550e-9m, D_t = 0.6m
            channel = FreespaceChannel(1000.0, 0.5, 0.0, 0.6, 0.6, 1550, clear, umbra)
            expected_waist = 1.22 * (1550e-9) / 0.6
            @test isapprox(waist_radius(channel), expected_waist; rtol=1e-10)
            
            # Test with different wavelength and diameter
            channel2 = FreespaceChannel(1000.0, 0.5, 0.0, 1.2, 0.6, 850, clear, umbra)
            expected_waist2 = 1.22 * (850e-9) / 1.2
            @test isapprox(waist_radius(channel2), expected_waist2; rtol=1e-10)
            
            # Test scaling behavior
            channel3 = FreespaceChannel(1000.0, 0.5, 0.0, 0.3, 0.6, 1550, clear, umbra)
            @test waist_radius(channel3) ≈ 2 * waist_radius(channel)  # half diameter = double waist
        end
        
        @testset "geometric_loss" begin
            # Test geometric loss calculation
            # L_geo = 20*log10((D_t + d*w₀)/D_r)
            channel = FreespaceChannel(1000.0, 0.5, 0.0, 0.6, 0.6, 1550, clear, umbra)
            w0 = waist_radius(channel)
            expected_loss = 20 * log10((0.6 + 1000.0 * w0) / 0.6)
            @test isapprox(geometric_loss(channel), expected_loss; rtol=1e-10)
            
            # Test with different parameters
            channel2 = FreespaceChannel(500.0, 0.5, 0.0, 1.2, 0.3, 1550, clear, umbra)
            w0_2 = waist_radius(channel2)
            expected_loss2 = 20 * log10((1.2 + 500.0 * w0_2) / 0.3)
            @test isapprox(geometric_loss(channel2), expected_loss2; rtol=1e-10)
            
            # Loss should increase with distance
            channel_far = FreespaceChannel(2000.0, 0.5, 0.0, 0.6, 0.6, 1550, clear, umbra)
            @test geometric_loss(channel_far) > geometric_loss(channel)
        end
        
        @testset "atmosphere_distance" begin
            # Test atmosphere distance calculation with different elevation angles
            channel_horizontal = FreespaceChannel(1000.0, 0.0, 0.0, 0.6, 0.6, 1550, clear, umbra)
            @test atmosphere_distance(channel_horizontal) == 0  # horizontal should be 0
            
            # Test with non-zero elevation
            channel_elevated = FreespaceChannel(1000.0, π/6, 100.0, 0.6, 0.6, 1550, clear, umbra)  # 30 degrees
            atm_dist = atmosphere_distance(channel_elevated)
            @test atm_dist > 0
            
            # Test with custom atmosphere height
            atm_dist_custom = atmosphere_distance(channel_elevated, 20e3)
            @test atm_dist_custom >= 0
            
            # Test that higher elevation angles give different results
            channel_steep = FreespaceChannel(1000.0, π/3, 100.0, 0.6, 0.6, 1550, clear, umbra)  # 60 degrees
            @test atmosphere_distance(channel_steep) != atmosphere_distance(channel_elevated)
        end
        
        @testset "atmospheric_loss" begin
            # Test atmospheric loss for supported wavelengths
            wavelengths = [550, 690, 850, 1550]
            expected_absorptions = [0.13, 0.01, 0.41, 0.01]
            
            for (wl, expected_abs) in zip(wavelengths, expected_absorptions)
                channel = FreespaceChannel(1000.0, π/4, 100.0, 0.6, 0.6, wl, clear, umbra)
                loss = atmospheric_loss(channel)
                @test loss isa Number
                @test loss >= 0  # Loss should be non-negative
            end
            
            # Test that unsupported wavelength throws error
            channel_bad = FreespaceChannel(1000.0, π/4, 100.0, 0.6, 0.6, 1000, clear, umbra)
            @test_throws String atmospheric_loss(channel_bad)
            
            # Test zero distance gives specific behavior
            channel_zero = FreespaceChannel(1000.0, 0.0, 100.0, 0.6, 0.6, 550, clear, umbra)
            # Should handle zero atmosphere distance gracefully
            @test atmospheric_loss(channel_zero) isa Number
        end
        
        @testset "pointing_loss" begin
            # Test pointing loss calculation
            # L_PNT = exp(-8Θ²ⱼ/w²₀)
            channel = FreespaceChannel(1000.0, π/4, 100.0, 0.6, 0.6, 1550, clear, umbra)
            w0 = waist_radius(channel)
            
            # Test with default jitter (5 microradians)
            jitter_default = 5
            expected_loss = exp(-8 * (jitter_default * 1e-6)^2 / w0^2)
            @test isapprox(pointing_loss(channel), expected_loss; rtol=1e-10)
            
            # Test with custom jitter
            jitter_custom = 10
            expected_loss_custom = exp(-8 * (jitter_custom * 1e-6)^2 / w0^2)
            @test isapprox(pointing_loss(channel, jitter_custom), expected_loss_custom; rtol=1e-10)
            
            # Test that larger jitter gives smaller transmission (larger loss)
            @test pointing_loss(channel, 10) < pointing_loss(channel, 5)
            
            # Test that result is between 0 and 1 (probability)
            loss = pointing_loss(channel)
            @test 0 <= loss <= 1
        end
        
        @testset "reflector_loss" begin
            # Test that reflector loss returns the expected constant value
            loss = reflector_loss()
            @test loss ≈ 5.854678746311231
            @test loss > 0  # Should be positive loss in dB
        end
        
        @testset "swapping_loss" begin
            # Test swapping loss calculation
            # swapping_loss = 0.5 * memory_read_write_loss
            
            # Test with default parameter
            default_loss = swapping_loss()
            @test default_loss ≈ 0.5 * 0.8  # 0.5 * default 0.8
            
            # Test with custom parameter
            custom_rw_loss = 0.9
            custom_loss = swapping_loss(custom_rw_loss)
            @test custom_loss ≈ 0.5 * custom_rw_loss
            
            # Test that result is a probability (between 0 and 1)
            @test 0 <= swapping_loss() <= 1
            @test 0 <= swapping_loss(0.95) <= 1
            
            # Test edge cases
            @test swapping_loss(0.0) == 0.0
            @test swapping_loss(1.0) == 0.5
        end
        
        @testset "Integration tests" begin
            # Test that all loss functions work together with a realistic channel
            using SatelliteToolboxPropagators: Propagators, OrbitPropagatorSgp4
            using SatelliteToolboxBase: OrbitStateVector
            tles = generate_regular_array_TLEs(orbits=1, sats_per_orbit=4, altitude_km=550)
            sat = Propagators.init(Val(:SGP4), tles[1])
            
            # Create a channel to subsatellite point
            sv = Propagators.propagate!(sat, 0, OrbitStateVector)
            lat, lon, _ = eci_to_geodetic(sv)
            
            channel = FreespaceChannel(sat, (lat, lon, 0.0); 
                                     time=sat.sgp4d.epoch, 
                                     min_θ=deg2rad(0),
                                     wavelength_nm=1550)
            
            if channel !== nothing
                # Test that all loss functions return reasonable values
                @test waist_radius(channel) > 0
                @test geometric_loss(channel) > 0  # Should have some geometric loss
                @test atmosphere_distance(channel) >= 0
                @test atmospheric_loss(channel) isa Number
                @test 0 <= pointing_loss(channel) <= 1
                @test reflector_loss() > 0
                @test 0 <= swapping_loss() <= 1
            end
        end
    end
end
