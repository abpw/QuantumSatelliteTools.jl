using QuantumSatelliteTools.Optimization
using QuantumSatelliteTools.GenerateGroundStations: generate_census_grid_gses
using QuantumSatelliteTools.GenerateSatellites: generate_regular_array_TLEs
using QuantumSatelliteTools.FreespaceChannels: FreespaceChannel, AbstractChannel
using QuantumSatelliteTools.LossCalculation: total_loss
using QuantumSatelliteTools.AstronomyGeometry: eci_to_ecef, ecef_to_geodetic,
    geodetic_to_ecef, seconds_per_day
using QuantumSatelliteTools.Simulator: Node, NodeLabel, Link, satellite, ground_station, simulate
using SatelliteToolboxBase: OrbitStateVector
using SatelliteToolboxPropagators: Propagators

struct TestLayout <: GroundStationLayout end
QuantumSatelliteTools.Optimization.generate_ground_stations(::TestLayout) = [(0.0, 0.0)]

struct TestChannel <: AbstractChannel
    transmissivity::Float64
end
QuantumSatelliteTools.LossCalculation.total_loss(channel::TestChannel) = 1 - channel.transmissivity

@testset "Optimization.jl" begin
    tle = only(generate_regular_array_TLEs(orbital_planes=1, sats_per_plane=1))
    sat = Node(Propagators.init(Val(:SGP4), tle), 0.0, NodeLabel(satellite, []))
    gses = [Node((0.0, lon), 0.0, NodeLabel(ground_station, [])) for lon in (0.0, 0.1, 0.2, 0.3)]
    channels = TestChannel.([0.5, 0.25, 0.125])
    links = [Link(sat, gs, channel) for (gs, channel) in zip(gses, channels)]

    @test DualDownlink(constrained=false)(links) == 0.5 * 0.25 + 0.5 * 0.125 + 0.25 * 0.125
    @test DualDownlink(constrained=true)(links) == 0.5 * 0.25
    dynamic_links = [Link(sat, gses[1], TestChannel(1.0)),
                     Link(sat, gses[2], TestChannel(1e-16))]
    @test DualDownlink(constrained=false)(dynamic_links) > 0
    tle2 = only(generate_regular_array_TLEs(name_prefix="OTHER", orbital_planes=1, sats_per_plane=1))
    sat2 = Node(Propagators.init(Val(:SGP4), tle2), 0.0, NodeLabel(satellite, []))
    links2 = vcat(links, [Link(sat2, gses[1], TestChannel(0.875)),
                         Link(gses[4], sat2, TestChannel(0.75))])
    @test DualDownlink(constrained=true)(links2) == 0.875 * 0.75 + 0.25 * 0.125

    @test generate_ground_stations(TestLayout()) == [(0.0, 0.0)]
    @test generate_census_grid_gses([(0.0, -0.001), (0.0, 0.0), (0.0, 0.001)]) == [(0.0, 0.0)]
    @test_throws ArgumentError generate_census_grid_gses([(0.0, 0.0)]; equatorial_distance_km=0)
    census = CensusGridLayout()
    stations = generate_ground_stations(census)
    @test !isempty(stations) && allunique(stations)
    empty!(stations)
    @test !isempty(generate_ground_stations(census))
    random_gses = generate_ground_stations(RandomLayout(3, 0))
    @test length(random_gses) == 3 && all(QuantumSatelliteTools.GenerateGroundStations.is_land, random_gses)
    @test random_gses == generate_ground_stations(RandomLayout(3, 0))
    @test length(QuantumSatelliteTools.Optimization.build_constellation(
        [Shell(2, 1, 10.0), Shell(3, 1, 20.0)], 550)) == 5
    @test simulate(Node[], 60, 30; evaluation_fn=_ -> 0, verbose=false) == [0, 0, 0]
    @test simulate(Node[], 0, 30; evaluation_fn=_ -> 0, verbose=false) == [0]

    future_time = sat.obj.sgp4d.epoch + 60 / seconds_per_day
    future_sat = Propagators.init(Val(:SGP4), tle)
    future_ecef = eci_to_ecef(Propagators.propagate!(future_sat, 60, OrbitStateVector).r; time=future_time)
    future_geo = ecef_to_geodetic(future_ecef)
    future_gs = (future_geo[1], future_geo[2])
    future_channel = FreespaceChannel(Propagators.init(Val(:SGP4), tle), future_gs;
                                      time=future_time, min_θ=0)
    expected_distance = sqrt(sum((geodetic_to_ecef(future_gs) .- future_ecef) .^ 2))
    @test future_channel.distance_m ≈ expected_distance

    shell = Shell(1, 1, 0.0)
    shell_sat = only(QuantumSatelliteTools.Optimization.build_constellation([shell], 550))
    shell_ecef = eci_to_ecef(Propagators.propagate!(shell_sat, 0, OrbitStateVector).r;
                             time=shell_sat.sgp4d.epoch)
    shell_geo = ecef_to_geodetic(shell_ecef)
    fixed_gses = [(shell_geo[1], mod(shell_geo[2] + offset + π, 2π) - π)
                  for offset in (-0.01, 0.0, 0.01)]
    direct = evaluate_constellation([shell], fixed_gses; duration_s=0, constrained=false,
                                    tx_aperture_m=0.1, rx_aperture_m=0.1)
    shell_node = Node(shell_sat, 0.0, NodeLabel(satellite, []))
    fixed_nodes = [Node(gs, 0.0, NodeLabel(ground_station, [])) for gs in fixed_gses]
    reference_links = [Link(shell_node, gs, FreespaceChannel(shell_sat, gs.obj;
                            time=shell_sat.sgp4d.epoch, transmitter_diameter_m=0.1,
                            receiver_diameter_m=0.1)) for gs in fixed_nodes]
    a, b, c = [1 - total_loss(link.channel) for link in reference_links]
    reference = a * b + a * c + b * c
    @test reference > 0
    @test direct ≈ reference * 1e9
    @test evaluate_constellation([shell], fixed_gses; duration_s=60, step_s=30,
        evaluation=_ -> 2.0, attempt_rate_hz=3) == 6
    for (duration_s, expected_samples) in ((seconds_per_day, 2880), (61, 3))
        samples = Ref(0)
        evaluate_constellation([shell], []; duration_s=duration_s,
            evaluation=_ -> (samples[] += 1))
        @test samples[] == expected_samples
    end
end
