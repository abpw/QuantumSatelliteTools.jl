using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumSatelliteTools.Simulator: build_optimized_constellation
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
    propagators = build_optimized_constellation()

    # Generate 100 population center ground stations.
    population_centers = generate_population_center_gses(100)

    # Generate ground station grids with 200, 400, and 800 km spacing.
    grid200 = generate_grid_gses(equatorial_distance_km=200)
    grid400 = generate_grid_gses(equatorial_distance_km=400)
    grid800 = generate_grid_gses(equatorial_distance_km=800)


end