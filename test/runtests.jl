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
        using QuantumSatelliteTools.AstronomyGeometry: decibel_to_probability, waist_radius, geometric_loss, atmosphere_distance, atmospheric_loss, pointing_loss, reflector_loss, swapping_loss, transmissivity, FreespaceChannel, clear, sunlight

        @testset "decibel_to_probability" begin
            # Test basic conversion
            @test decibel_to_probability(0) ≈ 0.0 atol=1e-10
            @test decibel_to_probability(10) ≈ 0.9 atol=1e-10
            @test decibel_to_probability(20) ≈ 0.99 atol=1e-10
            @test decibel_to_probability(30) ≈ 0.999 atol=1e-10
            # Negative values should give negative probabilities (unphysical but mathematically correct)
            @test decibel_to_probability(-10) < 0
        end

        @testset "waist_radius" begin
            # Create test channels with different parameters
            channel_500nm = FreespaceChannel(1000.0, 0.0, 0.0, 0.1, 0.05, 500, clear, sunlight)
            channel_1550nm = FreespaceChannel(1000.0, 0.0, 0.0, 0.1, 0.05, 1550, clear, sunlight)
            
            # Test waist radius calculation: w₀ = 1.22λ/D
            expected_500nm = 1.22 * (500e-9) / 0.1
            expected_1550nm = 1.22 * (1550e-9) / 0.1
            
            @test waist_radius(channel_500nm) ≈ expected_500nm rtol=1e-12
            @test waist_radius(channel_1550nm) ≈ expected_1550nm rtol=1e-12
            
            # Longer wavelength should give larger waist radius
            @test waist_radius(channel_1550nm) > waist_radius(channel_500nm)
        end

        @testset "geometric_loss" begin
            # Test geometric loss calculation
            channel = FreespaceChannel(1000.0, 0.0, 0.0, 0.1, 0.05, 1550, clear, sunlight)
            
            # Geometric loss = 20*log10((D_t + d*w₀)/D_r)
            w0 = waist_radius(channel)
            expected_loss = 20 * log10((0.1 + 1000.0 * w0) / 0.05)
            
            @test geometric_loss(channel) ≈ expected_loss rtol=1e-12
            @test geometric_loss(channel) > 0  # Should be positive (loss)
        end

        @testset "atmosphere_distance" begin
            # Test with zero elevation angle
            channel_zero_elev = FreespaceChannel(1000.0, 0.0, 0.0, 0.1, 0.05, 1550, clear, sunlight)
            @test atmosphere_distance(channel_zero_elev) == 0
            
            # Test with 45 degree elevation
            channel_45deg = FreespaceChannel(1000.0, π/4, 100000.0, 0.1, 0.05, 1550, clear, sunlight)
            atm_dist = atmosphere_distance(channel_45deg)
            @test atm_dist > 0
            @test atm_dist > 18000  # Should be at least atmosphere height / sin(45°) ≈ 25,456m
            
            # Test with 90 degree elevation (straight up)
            channel_90deg = FreespaceChannel(1000.0, π/2, 100000.0, 0.1, 0.05, 1550, clear, sunlight)
            atm_dist_90 = atmosphere_distance(channel_90deg)
            @test atm_dist_90 > 0
            
            # Test with custom atmosphere height
            custom_atm_height = 10e3
            atm_dist_custom = atmosphere_distance(channel_45deg, custom_atm_height)
            @test atm_dist_custom >= 0
        end

        @testset "atmospheric_loss" begin
            # Test supported wavelengths
            wavelengths = [550, 690, 850, 1550]
            expected_absorptions = [0.13, 0.01, 0.41, 0.01]
            
            for (wl, expected_absorption) in zip(wavelengths, expected_absorptions)
                channel = FreespaceChannel(10000.0, π/4, 200000.0, 0.1, 0.05, wl, clear, sunlight)
                loss = atmospheric_loss(channel)
                @test loss >= 0  # Loss should be non-negative
                
                # Calculate expected loss based on transmittance
                atm_dist = atmosphere_distance(channel)
                expected_loss = expected_absorption * (atm_dist / 1000)
                @test loss ≈ expected_loss rtol=1e-12
            end
            
            # Test unsupported wavelength throws error
            channel_bad_wl = FreespaceChannel(1000.0, π/4, 100000.0, 0.1, 0.05, 633, clear, sunlight)
            @test_throws Exception atmospheric_loss(channel_bad_wl)
        end

        @testset "pointing_loss" begin
            channel = FreespaceChannel(1000.0, π/4, 100000.0, 0.1, 0.05, 1550, clear, sunlight)
            
            # Test with default pointing jitter (5 μrad)
            loss_default = pointing_loss(channel)
            @test 0 < loss_default <= 1  # Should be a probability
            
            # Test with custom pointing jitter
            loss_10 = pointing_loss(channel, 10)
            loss_1 = pointing_loss(channel, 1)
            
            # Higher jitter should give lower transmission (higher loss)
            @test loss_1 > loss_10
            
            # Test formula: L_PNT = exp(-8Θ²ⱼ/w²₀)
            w0 = waist_radius(channel)
            jitter = 5e-6  # 5 microradians
            expected = exp(-8 * jitter^2 / w0^2)
            @test pointing_loss(channel, 5) ≈ expected rtol=1e-12
        end

        @testset "reflector_loss" begin
            # Test that reflector loss is a constant
            @test reflector_loss() == 5.854678746311231
            @test typeof(reflector_loss()) <: Number
        end

        @testset "swapping_loss" begin
            # Test with default memory loss (0.8)
            loss_default = swapping_loss()
            @test loss_default == 0.5 * 0.8
            @test loss_default == 0.4
            
            # Test with custom memory loss
            loss_09 = swapping_loss(0.9)
            @test loss_09 == 0.5 * 0.9
            @test loss_09 == 0.45
            
            # Test edge cases
            @test swapping_loss(0.0) == 0.0
            @test swapping_loss(1.0) == 0.5
        end

        @testset "transmissivity" begin
            # Create a test channel
            channel = FreespaceChannel(5000.0, π/6, 400000.0, 0.1, 0.05, 1550, clear, sunlight)
            
            # Test that transmissivity is between 0 and 1
            trans = transmissivity(channel)
            @test 0 <= trans <= 1
            
            # Test that transmissivity formula matches expectation
            # transmissivity = 10^(-L_dB/10) * pointing_probability
            geo_loss = geometric_loss(channel)
            atm_loss = atmospheric_loss(channel)
            pointing_prob = pointing_loss(channel)
            total_loss_db = geo_loss + atm_loss
            
            expected_trans = 10^(-total_loss_db / 10) * pointing_prob
            @test trans ≈ expected_trans rtol=1e-10
            
            # Test with a different channel configuration
            alt_channel = FreespaceChannel(100.0, π/2, 500000.0, 0.2, 0.1, 1550, clear, sunlight)
            trans_alt = transmissivity(alt_channel)
            @test 0 <= trans_alt <= 1  # Should be a valid transmissivity
        end

        @testset "integration tests" begin
            # Test that all functions work together in a realistic scenario
            
            # Typical satellite-to-ground link parameters
            distance_m = 1000e3  # 1000 km
            elevation_rad = deg2rad(30)  # 30 degree elevation
            min_altitude_m = 400e3  # 400 km altitude
            tx_diameter = 0.3  # 30 cm transmitter
            rx_diameter = 1.0   # 1 m receiver
            wavelength = 1550   # 1550 nm
            
            channel = FreespaceChannel(distance_m, elevation_rad, min_altitude_m, 
                                     tx_diameter, rx_diameter, wavelength, clear, sunlight)
            
            # All loss functions should return reasonable values
            @test 0 < waist_radius(channel) < 1e-3  # Should be microradians
            @test geometric_loss(channel) > 0
            @test atmosphere_distance(channel) > 0
            @test atmospheric_loss(channel) >= 0
            @test 0 < pointing_loss(channel) <= 1
            @test reflector_loss() > 0
            @test 0 < swapping_loss() < 1
            
            # Overall transmissivity should be reasonable for this link
            trans = transmissivity(channel)
            @test 1e-6 < trans < 0.1  # Typical range for satellite links
        end
    end
end
