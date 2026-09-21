using Dates: DateTime
using Test: @testset, @test

using SatelliteToolboxTle: read_tle
using SatelliteToolboxPropagators: Propagators

using QuantumSatelliteTools
using QuantumSatelliteTools.constants
using QuantumSatelliteTools.types
using QuantumSatelliteTools.FreespaceChannels

@testset "FreespaceChannel" begin
    @testset "direct construction with defaults" begin
        distance = 1_000_000  # 1000 km
        elevation = deg2rad(45)

        ch = FreespaceChannel(distance, elevation)
        @test ch isa FreespaceChannel
        @test ch.distance_m == distance
        @test ch.min_altitude_m == 0
        @test ch.elevation_angle_rad == elevation
        @test ch.transmitter_diameter_m == 0.6
        @test ch.receiver_diameter_m == 0.6
        @test ch.wavelength_nm == 1550
        @test ch.conditions == clear
        @test ch.light_condition == umbra
    end

    @testset "direct construction with keywords" begin
        distance = 500000.0
        elevation = deg2rad(30.0)

        ch = FreespaceChannel(
            distance,
            elevation;
            min_altitude_m=100,
            transmitter_diameter_m=0.5,
            receiver_diameter_m=0.8,
            wavelength_nm=850,
            conditions=rain,
            light_condition=penumbra,
        )
        @test ch isa FreespaceChannel
        @test ch.distance_m == distance
        @test ch.elevation_angle_rad == elevation
        @test ch.min_altitude_m == 100
        @test ch.transmitter_diameter_m == 0.5
        @test ch.receiver_diameter_m == 0.8
        @test ch.wavelength_nm == 850
        @test ch.conditions == rain
        @test ch.light_condition == penumbra
    end

    @testset "sat-gs construction" begin
        tle_string = TEST_TLES[:ISS]

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)

        ch = FreespaceChannel(sat, gs)

        if ch !== nothing
            @test ch isa FreespaceChannel

            @test 6.7e6 < ch.distance_m < 2.5e7 # ISS distance ∈ [6_700, 25_000] km
            @test 0 <= ch.elevation_angle_rad <= π/2
        end
    end

    @testset "sat-gs with min_angle" begin
        tle_string = TEST_TLES[:ISS]

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)

        # Low minimum elevation angle - likely visible
        ch_low_θ = FreespaceChannel(sat, gs; time=sat.sgp4d.epoch, min_angle=deg2rad(10))
        @test ch_low_θ === nothing || ch_low_θ isa FreespaceChannel
        if ch_low_θ !== nothing
            @test ch_low_θ.elevation_angle_rad >= deg2rad(10)
        end

        # High minimum elevation angle - may not be visible
        ch_high_θ = FreespaceChannel(sat, gs; time=sat.sgp4d.epoch, min_angle=deg2rad(60))
        @test ch_high_θ === nothing || ch_high_θ isa FreespaceChannel
        if ch_high_θ !== nothing
            @test ch_high_θ.elevation_angle_rad >= deg2rad(60)
        end
    end

    @testset "sat-gs construction with keywords" begin
        tle_string = TEST_TLES[:ISS]

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)

        ch = FreespaceChannel(
            sat,
            gs;
            time=DateTime(2026, 8, 13, 12, 0, 0),
            transmitter_diameter_m=0.5,
            receiver_diameter_m=0.8,
            wavelength_nm=850,
            conditions=rain,
            min_angle=deg2rad(10),
        )

        if ch !== nothing
            @test ch isa FreespaceChannel
            @test ch.min_altitude_m == 100
            @test ch.transmitter_diameter_m == 0.5
            @test ch.receiver_diameter_m == 0.8
            @test ch.wavelength_nm == 850
            @test ch.conditions == rain
        end
    end

    @testset "sat-sat construction" begin
        tle1_string = TEST_TLES[:ISS]
        tle2_string = TEST_TLES[:Hubble]

        tle1 = read_tle(tle1_string, verify_checksum=false)
        tle2 = read_tle(tle2_string, verify_checksum=false)

        sat1 = Propagators.init(Val(:SGP4), tle1)
        sat2 = Propagators.init(Val(:SGP4), tle2)

        ch = FreespaceChannel(sat1, sat2; time=sat1.sgp4d.epoch)

        # May be nothing if obstructed or not visible
        if ch !== nothing
            @test ch isa FreespaceChannel

            # ISS altitude ~400 km, Hubble altitude ~560 km
            # Minimum distance: ~160 km (if directly above/below)
            # Maximum distance: ~20,000 km (antipodal on orbits)
            @test 1.6e5 < ch.distance_m < 2.5e7  # Distance ∈ [160, 25_000] km
        end
    end

    @testset "sat-sat construction with keywords" begin
        tle1_string = TEST_TLES[:ISS]
        tle2_string = TEST_TLES[:Hubble]

        tle1 = read_tle(tle1_string, verify_checksum=false)
        tle2 = read_tle(tle2_string, verify_checksum=false)

        sat1 = Propagators.init(Val(:SGP4), tle1)
        sat2 = Propagators.init(Val(:SGP4), tle2)

        ch = FreespaceChannel(
            sat1,
            sat2;
            time=DateTime(2026, 8, 13, 12, 0, 0),
            transmitter_diameter_m=0.5,
            receiver_diameter_m=0.8,
            wavelength_nm=850,
            conditions=rain,
        )

        if ch !== nothing
            @test ch isa FreespaceChannel
            @test ch.min_altitude_m == 100
            @test ch.transmitter_diameter_m == 0.5
            @test ch.receiver_diameter_m == 0.8
            @test ch.wavelength_nm == 850
            @test ch.conditions == rain
        end
    end

    @testset "atmospheric conditions" begin
        distance = 500000.0
        elevation = deg2rad(45.0)

        conditions = [clear, rain, fog, snow]

        channels = [
            FreespaceChannel(distance, elevation; conditions=cond) for cond in conditions
        ]

        @test all(ch isa FreespaceChannel for ch in channels)
        @test all(ch.conditions == cond for (ch, cond) in zip(channels, conditions))
    end

    @testset "datetime input" begin
        tle_string = TEST_TLES[:ISS]

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)
        dt = DateTime(2026, 8, 13, 12, 0, 0)

        ch = FreespaceChannel(sat, gs; time=dt)

        @test ch === nothing || ch isa FreespaceChannel
    end
end