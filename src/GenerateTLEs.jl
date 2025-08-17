using SatelliteToolboxTle
using SatelliteAnalysis
using Downloads
using CSV
using ..AstronomyGeometry: semimajor_radius, earth_mass_kg, seconds_per_day, G

const satcat_file_path = joinpath(@__DIR__, "../databases/satcat.csv")
const starlink_TLE_file_path = joinpath(@__DIR__, "../databases/starlink.tle")

"""
Download the latest `satcat.csv` catalog from Celestrak and save it locally.
"""
function refresh_satcat()
    Downloads.download("https://celestrak.org/pub/satcat.csv", satcat_file_path)
end

"""
Update the stored Starlink TLEs by querying Celestrak for new entries.
    Minimizes calls to Celestrak by only querying for Starlink satellites
    that are listed in satcat.csv as active but do not appear in the 
    local Starlink TLE database.

TODO: remove inactive satellites

# Keyword Arguments
- `verbose::Bool=false`: If true, prints status messages during fetching.
"""
function update_starlink_tles(;verbose::Bool=false)
    existing_IDs = []
    if isfile(starlink_TLE_file_path)
        existing_TLEs = SatelliteToolboxTle.read_tles_from_file(starlink_TLE_file_path)
        existing_IDs = [tle.satellite_number for tle in existing_TLEs]
    end
    save_TLEs(get_active_satellite_TLEs(OBJECT_NAME="starlink", row_filter=row->!(parse(Int, row.NORAD_CAT_ID)∈existing_IDs), verbose=verbose), starlink_TLE_file_path)
end

"""
Save a vector of TLE objects to a specified file.

# Arguments
- `TLE_vector::Vector`: The vector of TLE objects to save (see [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)).
- `file_path::String`: The path to the output file.
"""
function save_TLEs(TLE_vector::Vector, file_path::String)
    TLEs_str = ""
    for TLE_obj in TLE_vector
        TLEs_str *= convert(String, TLE_obj) * "\n"
    end
    open(file->write(file, TLEs_str), file_path, create=true, write=true)
end

const tle_fetcher = create_tle_fetcher(CelestrakTleFetcher)

"""
Retrieve TLEs of active payload satellites from `satcat.csv` and Celestrak.

# Keyword Arguments
- `OBJECT_NAME::String=""`: Optional filter to match part of the satellite name.
- `row_filter::Function=sat_row->true`: Function to filter rows from `satcat.csv`.
- `return_partial::Bool=true`: If true, returns partial results even if some fetches fail.
- `verbose=false`: If true, prints progress info.

# Returns
- `Vector{TLE}`: A vector of fetched TLE objects.
"""
function get_active_satellite_TLEs(;OBJECT_NAME::String="", row_filter::Function=sat_row->true, return_partial::Bool=true, verbose=false)
    NORAD_CAT_IDs = []
    OBJECT_NAME = uppercase(OBJECT_NAME)
    for row ∈ CSV.Rows(satcat_file_path)
        if typeof(row.OPS_STATUS_CODE) != Missing && only(row.OPS_STATUS_CODE) == '+' && row.OBJECT_TYPE == "PAY" && occursin(OBJECT_NAME, row.OBJECT_NAME) && row_filter(row)
            push!(NORAD_CAT_IDs, parse(Int, row.NORAD_CAT_ID))
        end
    end
    TLEs = []
    for id_number in NORAD_CAT_IDs
        try
            push!(TLEs, fetch_tles(tle_fetcher, satellite_number=id_number)[1])
        catch
            if verbose
                if length(TLEs) > 0
                    println("Fetched $(length(TLEs)) TLEs before being denied. Last fetched ID#: $(TLEs[end].satellite_number)")
                else
                    println("Denied by Celestrak before fetching any TLEs")
                end
            end
            if return_partial
                return TLEs
            end
            return []
        end
    end
    return TLEs
end

"""
Generate synthetic TLEs for a regular constellation of satellites.

# Keyword Arguments
- `name_prefix::String="SAT"`: Prefix for satellite names.
- `orbits::Int=12`: Number of orbital planes.
- `sats_per_orbit::Int=18`: Satellites per orbital plane.
- `altitude_km::Int=500`: Altitude of the orbits in kilometers.
- `inclination_rad::Float64=π/2`: Inclination angle in radians.
- `frozen_orbits::Bool=false`: Whether to use frozen orbit parameters. See https://juliaspace.github.io/SatelliteAnalysis.jl/stable/man/frozen_orbits/.

# Returns
- `Vector{TLE}`: Vector of synthetic TLEs.
"""
function generate_regular_array_TLEs(;
            name_prefix::String="SAT",
            orbits::Int=12, sats_per_orbit::Int=18,
            altitude_km::Int=500,
            inclination_rad::Float64=π/2,
            frozen_orbits::Bool=false)

    TLEs = []
    for orbit ∈ 1:orbits
        Omega = 180*(orbit - 1)/orbits
        for sat_in_orbit ∈ 1:sats_per_orbit
            anomaly = 360/sats_per_orbit*(sat_in_orbit - 1 + (orbit - 1)/orbits)
            number = (orbit - 1)*sats_per_orbit + sat_in_orbit
            name = name_prefix * " " * string(number)
            if frozen_orbits
                eccentricity, argument_of_perigee = SatelliteAnalysis.frozen_orbit(altitude_km*1000+semimajor_radius, rad2deg(inclination_rad))
            else
                eccentricity, argument_of_perigee = (10^-7, 90)
            end
            push!(TLEs, TLE(name=name,
                            epoch_year=25,
                            epoch_day=1,
                            inclination=inclination_rad/π*180,
                            raan=Omega,
                            eccentricity=eccentricity,
                            argument_of_perigee=argument_of_perigee,
                            mean_anomaly=anomaly,
                            mean_motion=seconds_per_day*√(G*earth_mass_kg/(altitude_km*1000+semimajor_radius)^3)/2/π))
        end
    end
    return TLEs
end