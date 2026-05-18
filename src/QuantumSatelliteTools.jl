module QuantumSatelliteTools

using GeoDatasets
using SatelliteToolboxTle, SatelliteToolboxPropagators, SatelliteToolboxTransformations, SatelliteToolboxBase
using SatelliteAnalysis
using Downloads, CSV, Random, Dates, StaticArrays

export AstronomyGeometry, GenerateSatellites, GenerateGroundStations, FreespaceChannels, LossCalculation

module AstronomyGeometry
export semimajor_radius, semiminor_radius, G, earth_mass_kg, seconds_per_day, equatorial_circumference_km, sin_60
export GS, norm, ellipsoid_line_intersection
export geodetic_to_ecef, ecef_to_geodetic, eci_to_ecef, eci_to_geodetic
export gs_gs_distance, sat_sat_distance
export LightCondition, sunlight, penumbra, umbra
export sun_position, light_at_point
export Conditions, clear, fog, rain, snow
include("AstronomyGeometry.jl")
end

module GenerateGroundStations
export is_within_min_distance, is_land
export generate_population_center_gses, generate_city_gses, generate_random_gses, generate_equispaced_gses, generate_grid_gses
include("GenerateGroundStations.jl")
end

module GenerateSatellites
export refresh_satcat, update_starlink_tles, save_TLEs, get_active_satellite_TLEs, generate_regular_array_TLEs
include("GenerateSatellites.jl")
end

module FreespaceChannels
export FreespaceChannel
include("FreespaceChannels.jl")
end

module LossCalculation
export atmospheric_loss_dB, geometric_loss_dB, pointing_loss_dB, reflector_loss, swapping_loss, total_loss
include("LossCalculation.jl")
end

module Simulator
export build_optimized_constellation, simulate_step
include("Simulator.jl")
end

end
