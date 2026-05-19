using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumSatelliteTools.Simulator: build_optimized_constellation, Node, NodeLabels, simulate
using QuantumSatelliteTools.GenerateGroundStations: generate_grid_gses, generate_population_center_gses

"""
Perform the grid ground station selection experiment.

# Arguments
- `duration`: total simulation duration in seconds.
- `step`: time step in seconds.

# Returns
- Nothing
"""
function main(;
    duration_s::Int = 30,
    step_s::Int = 30,
)
    # Build a constellation of 1000 satellites based on optimized parameters
    satellites = build_optimized_constellation()

    # Generate 100 population center ground stations.
    population_center_gses = generate_population_center_gses(100)

    # Generate ground station grids with 200, 400, and 800 km spacing.
    grid_gses = Dict(
        200 => generate_grid_gses(equatorial_distance_km=200),
        400 => generate_grid_gses(equatorial_distance_km=400),
        800 => generate_grid_gses(equatorial_distance_km=800),
    )

    # perform experiments for each grid scenario
    for (spacing, aux_gses) in grid_gses
        println("Simulating scenario with grid spacing of $spacing km...")

        # Construct the nodes for this scenario
        groups = [
            (satellites, _ -> 0.0, [NodeLabels.sat]),
            (population_center_gses, _ -> 1.0, [NodeLabels.gs]),
            (aux_gses, _ -> 0.0, [NodeLabels.aux_gs]),
        ]
        nodes = [
            Node(obj, weight_fn(obj), labels)
            for (objs, weight_fn, labels) in groups
            for obj in objs
        ]

        # Simulate the scenario and collect metrics
        metrics = simulate(nodes, duration_s, step_s)

        println("Metrics for grid spacing $spacing km: $metrics")
    end

end