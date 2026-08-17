using Dates: DateTime
using Test: @testset, @test

using SatelliteToolboxTle: read_tle
using SatelliteToolboxPropagators: Propagators

using QuantumSatelliteTools
using QuantumSatelliteTools.FreespaceChannels

@testset "FreespaceChannel" begin
    @testset "direct construction" begin
        distance = 1000000.0  # 1000 km
        elevation = deg2rad(45.0)

        ch = FreespaceChannel(distance, elevation)
        @test ch isa FreespaceChannel
        @test ch.distance_m == distance
        @test ch.elevation_angle_rad == elevation
    end

    @testset "direct construction with keywords" begin
        distance = 500000.0
        elevation = deg2rad(30.0)

        ch = FreespaceChannel(
            distance,
            elevation;
            transmitter_diameter_m=0.5,
            receiver_diameter_m=0.8,
            wavelength_nm=1064,
        )
        @test ch isa FreespaceChannel
        @test ch.distance_m == distance
    end

    @testset "sat-gs construction" begin
        tle_string = """ISS (ZARYA)
            1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
            2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)  # Equator, prime meridian, sea level

        ch = FreespaceChannel(sat, gs; time=sat.sgp4d.epoch)

        # May be nothing if not visible, or a FreespaceChannel
        if ch !== nothing
            @test ch isa FreespaceChannel
            @test ch.distance_m > 0
            @test ch.elevation_angle_rad >= -π/2
        end
    end

    @testset "sat-gs with min_angle" begin
        tle_string = """ISS (ZARYA)
            1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
            2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)

        # High minimum elevation angle - may not be visible
        ch = FreespaceChannel(sat, gs; time=sat.sgp4d.epoch, min_angle=deg2rad(60))

        @test ch === nothing || ch isa FreespaceChannel
    end

    @testset "sat-sat construction" begin
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

        ch = FreespaceChannel(sat1, sat2; time=sat1.sgp4d.epoch)

        # May be nothing if obstructed or not visible
        if ch !== nothing
            @test ch isa FreespaceChannel
            @test ch.distance_m > 0
        end
    end

    @testset "return type or nothing" begin
        distance = 1000000.0
        elevation = deg2rad(30.0)

        ch = FreespaceChannel(distance, elevation)
        @test ch isa FreespaceChannel
    end

    @testset "atmospheric conditions" begin
        distance = 500000.0
        elevation = deg2rad(45.0)

        ch_clear = FreespaceChannel(distance, elevation; conditions=clear)
        ch_rain = FreespaceChannel(distance, elevation; conditions=rain)
        ch_fog = FreespaceChannel(distance, elevation; conditions=fog)
        ch_snow = FreespaceChannel(distance, elevation; conditions=snow)

        @test ch_clear isa FreespaceChannel
        @test ch_rain isa FreespaceChannel
        @test ch_fog isa FreespaceChannel
        @test ch_snow isa FreespaceChannel
    end

    @testset "datetime input" begin
        tle_string = """ISS (ZARYA)
            1 25544U 98067A   26225.00000000  .00016717  00000-0  29286-3 0  9990
            2 25544  51.6412 339.8014 0004346  86.3973 273.7418 15.54115342431000"""

        tle = read_tle(tle_string, verify_checksum=false)
        sat = Propagators.init(Val(:SGP4), tle)
        gs = (0.0, 0.0, 0.0)
        dt = DateTime(2026, 8, 13, 12, 0, 0)

        ch = FreespaceChannel(sat, gs; time=dt)

        @test ch === nothing || ch isa FreespaceChannel
    end

    @testset "different telescopes" begin
        distance = 500000.0
        elevation = deg2rad(45.0)

        ch1 = FreespaceChannel(distance, elevation; transmitter_diameter_m=0.1)
        ch2 = FreespaceChannel(distance, elevation; transmitter_diameter_m=1.0)

        @test ch1 isa FreespaceChannel
        @test ch2 isa FreespaceChannel
    end

    @testset "different wavelengths" begin
        distance = 500000.0
        elevation = deg2rad(45.0)

        ch_ir = FreespaceChannel(distance, elevation; wavelength_nm=1550)
        ch_vis = FreespaceChannel(distance, elevation; wavelength_nm=532)

        @test ch_ir isa FreespaceChannel
        @test ch_vis isa FreespaceChannel
    end
end