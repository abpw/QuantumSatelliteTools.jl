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

export constants
export types
export helpers

export AstronomyGeometry
export FreespaceChannels
export GenerateGroundStations
export GenerateSatellites
export LossCalculation
export Simulator
export VisualizeMap

module constants
export AU_M
export SEMIMAJOR_RADIUS, SEMIMINOR_RADIUS
export EQUATORIAL_CIRCUMFERENCE_KM
export EARTH_MASS_KG
export G
export SECONDS_PER_DAY
export SIN_60
export LAND_SEA_MASK
include("constants.jl")
end

module types
export InputGS, GS, Point3D
export AbstractChannel
export Conditions, clear, fog, rain, snow
export LightCondition, umbra, penumbra, sunlight
export Propagator
export NodeTypes, satellite, ground_station
export NodeLabel, Node, Link
export Time
include("types.jl")
end

module helpers
export norm
export dB_to_prob, prob_to_dB
export validate_non_negative_finite, validate_coordinates
export is_satellite, is_ground_station
export is_duplicate, is_land
include("helpers.jl")
end

module AstronomyGeometry
export ellipsoid_line_intersection
export geodetic_to_ecef, ecef_to_geodetic, eci_to_ecef, eci_to_geodetic
export gs_gs_distance, sat_sat_distance
export sun_position, light_at_point
include("AstronomyGeometry.jl")
end

module FreespaceChannels
export FreespaceChannel
include("FreespaceChannels.jl")
end

module GenerateGroundStations
export is_outside_min_distance
export generate_population_center_gses
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
export atmospheric_loss_dB, geometric_loss_dB, pointing_loss_dB
export reflector_loss
export swapping_loss, swapping_loss_dB
export total_loss, total_loss_dB
include("LossCalculation.jl")
end

module Simulator
export build_optimized_constellation
export simulate
include("Simulator.jl")
end

module VisualizeMap
export plot_gses, plot_gs_selections
include("VisualizeMap.jl")
end
end
