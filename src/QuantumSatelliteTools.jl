module QuantumSatelliteTools

using CSV
using Dates
using Downloads
using Random
using StaticArrays

using GeoDatasets
using SatelliteAnalysis
using SatelliteToolboxBase
using SatelliteToolboxPropagators
using SatelliteToolboxTle
using SatelliteToolboxTransformations

export AstronomyGeometry
export FreespaceChannels
export GenerateGroundStations
export GenerateSatellites
export LossCalculation
export Simulator
export VisualizeMap

module AstronomyGeometry
export semimajor_radius, semiminor_radius, G, earth_mass_kg, equatorial_circumference_km
export seconds_per_day, sin_60
export GS, norm, ellipsoid_line_intersection
export geodetic_to_ecef, ecef_to_geodetic, eci_to_ecef, eci_to_geodetic
export gs_gs_distance, sat_sat_distance
export LightCondition, sunlight, penumbra, umbra
export sun_position, light_at_point
export Conditions, clear, fog, rain, snow
include("AstronomyGeometry.jl")
end

module FreespaceChannels
export FreespaceChannel
include("FreespaceChannels.jl")
end

module GenerateGroundStations
export is_within_min_distance, is_land
export generate_population_center_gses, _generate_population_center_gses
export generate_city_gses
export generate_random_gses
export generate_equispaced_gses
export generate_grid_gses
include("GenerateGroundStations.jl")
end

module GenerateSatellites
export refresh_satcat, update_starlink_TLEs, save_TLEs
export get_active_satellites, get_active_satellite_TLEs
export generate_regular_array, generate_regular_array_TLEs
include("GenerateSatellites.jl")
end

module LossCalculation
export dB_to_prob, prob_to_dB
export atmospheric_loss_dB, geometric_loss_dB, pointing_loss_dB
export reflector_loss
export swapping_loss, swapping_loss_dB
export total_loss, total_loss_dB
include("LossCalculation.jl")
end

module Simulator
export build_optimized_constellation
export Node, Link, NodeLabel
export GroundStation, Propagator
export satellite, ground_station
export simulate
include("Simulator.jl")
end

module VisualizeMap
export plot_gses, plot_gs_selections
include("VisualizeMap.jl")
end

end
