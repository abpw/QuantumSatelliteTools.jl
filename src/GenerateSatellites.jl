using Downloads: download

using CSV: Rows
using SatelliteAnalysis: frozen_orbit
using SatelliteToolboxPropagators: Propagators
using SatelliteToolboxTle: TLE, CelestrakTleFetcher, create_tle_fetcher, read_tles_from_file

using ..constants: SEMIMAJOR_RADIUS, EARTH_MASS_KG, SECONDS_PER_DAY, G

const SATCAT_FILE_PATH = joinpath(@__DIR__, "../databases/satcat.csv")
const STARLINK_TLE_FILE_PATH = joinpath(@__DIR__, "../databases/starlink.tle")

"""
    refresh_satcat()

Download the latest `satcat.csv` catalog from Celestrak and save it locally.
"""
function refresh_satcat()
    download("https://celestrak.org/pub/satcat.csv", SATCAT_FILE_PATH)
end

"""
    update_starlink_TLEs(; verbose::Bool=false)

Update stored Starlink TLEs by querying Celestrak for new entries.

Minimizes calls to Celestrak by only querying for Starlink satellites that are listed in
satcat.csv as active but do not appear in the local Starlink TLE database.

# Keyword Arguments

  - `verbose::Bool=false`: If `true`, prints status messages during fetching.

# Notes

  - TODO: remove inactive satellites
"""
function update_starlink_TLEs(; verbose::Bool=false)
    existing_IDs::Vector{Int} = []
    if isfile(STARLINK_TLE_FILE_PATH)
        existing_TLEs = read_tles_from_file(STARLINK_TLE_FILE_PATH)
        if verbose
            println("Using $(length(existing_TLEs)) existing TLEs")
        end
        existing_IDs = [tle.satellite_number for tle in existing_TLEs]
    end
    save_TLEs(
        get_active_satellite_TLEs(;
            OBJECT_NAME="starlink",
            row_filter=row->!(parse(Int, row.NORAD_CAT_ID)∈existing_IDs),
            verbose=verbose,
        ),
        STARLINK_TLE_FILE_PATH,
    )
end

"""
    save_TLEs(TLE_vector::Vector{TLE}, file_path::String) -> Nothing

Save vector of TLE objects to a specified file.

# Arguments

  - `TLE_vector::Vector`: TLE objects to save
    (see [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl)).
  - `file_path::String`: Path to the output file.
"""
function save_TLEs(TLE_vector::Vector{TLE}, file_path::String)
    TLEs_str = join(string.(TLE_vector), "\n")
    open(file->write(file, TLEs_str), file_path; create=true, write=true)
end

const TLE_FETCHER = create_tle_fetcher(CelestrakTleFetcher)

"""
    get_active_satellite_TLEs(; <keyword arguments>)

Retrieve TLEs of active payload satellites from `satcat.csv` and Celestrak.

# Keyword Arguments

  - `OBJECT_NAME::String=""`: Optional filter to match part of the satellite name.
  - `row_filter::Function=sat_row->true`: Function to filter rows from `satcat.csv`.
  - `return_partial::Bool=true`: If `true`, returns partial results even if some fetches
    fail.
  - `verbose::Bool=false`: If `true`, prints progress info.

# Returns

  - `Vector{TLE}`: Fetched TLE objects.
"""
function get_active_satellite_TLEs(;
    OBJECT_NAME::String="",
    row_filter::Function=sat_row->true,
    return_partial::Bool=true,
    verbose::Bool=false,
)
    NORAD_CAT_IDs::Vector{Int} = []
    OBJECT_NAME = uppercase(OBJECT_NAME)
    for row in CSV.Rows(SATCAT_FILE_PATH)
        if !ismissing(row.OPS_STATUS_CODE) &&
            only(row.OPS_STATUS_CODE) == '+' &&
            row.OBJECT_TYPE == "PAY" &&
            occursin(OBJECT_NAME, row.OBJECT_NAME) &&
            row_filter(row)
            push!(NORAD_CAT_IDs, parse(Int, row.NORAD_CAT_ID))
        end
    end
    TLEs = []
    for id_number in NORAD_CAT_IDs
        try
            push!(TLEs, fetch_tles(TLE_FETCHER; satellite_number=id_number)[1])
        catch
            if verbose
                if length(TLEs) > 0
                    println(
                        "Fetched $(length(TLEs)) TLEs before being denied. " *
                            "Last fetched ID#: $(TLEs[end].satellite_number)",
                    )
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
    get_active_satellites(; <keyword arguments>)

Retrieve and return Propagator objects of active payload satellites from `satcat.csv` and
Celestrak.

# Keyword Arguments

  - `OBJECT_NAME::String=""`: Optional filter to match part of the satellite name.
  - `row_filter::Function=sat_row->true`: Function to filter rows from `satcat.csv`.
  - `return_partial::Bool=true`: If `true`, returns partial results even if some fetches
    fail.
  - `verbose::Bool=false`: If `true`, prints progress info.

# Returns

  - `Vector{Propagator}`: Fetched Propagator objects.
"""
function get_active_satellites(;
    OBJECT_NAME::String="",
    row_filter::Function=sat_row->true,
    return_partial::Bool=true,
    verbose::Bool=false,
)
    TLEs = get_active_satellite_TLEs(;
        OBJECT_NAME=OBJECT_NAME,
        row_filter=row_filter,
        return_partial=return_partial,
        verbose=verbose,
    )
    return [Propagators.init(Val(:SGP4), tle) for tle in TLEs]
end

"""
    generate_regular_array_TLEs(; <keyword arguments>)

Generate synthetic TLEs for a regular constellation of satellites.

# Keyword Arguments

  - `name_prefix::String="SAT"`: Prefix for satellite names.
  - `orbital_planes::Int=12`: Number of orbital planes.
  - `sats_per_plane::Int=18`: Number of satellites per orbital plane.
  - `altitude_km::Int=500`: Altitude of the orbits in kilometers.
  - `inclination_rad::Float64=π/2`: Inclination angle in radians.
  - `frozen_orbits::Bool=false`: Whether to use frozen orbit parameters (see
    [SatelliteAnalysis.jl/frozen_orbits](https://juliaspace.github.io/SatelliteAnalysis.jl/stable/man/frozen_orbits/)).

# Returns

  - `Vector{TLE}`: Synthetic TLEs.
"""
function generate_regular_array_TLEs(;
    name_prefix::String="SAT",
    orbital_planes::Int=12,
    sats_per_plane::Int=18,
    altitude_km::Int=500,
    inclination_rad::Float64=π/2,
    frozen_orbits::Bool=false,
)
    TLEs::Vector{TLE} = []
    for plane in 1:orbital_planes
        Ω = 180 * (plane - 1)/orbital_planes
        for sat_in_orbit in 1:sats_per_plane
            anomaly = 360/sats_per_plane*(sat_in_orbit - 1 + (plane - 1)/orbital_planes)
            number = (plane - 1)*sats_per_plane + sat_in_orbit
            name = name_prefix * " " * string(number)
            if frozen_orbits
                eccentricity, argument_of_perigee = frozen_orbit(
                    altitude_km*1000+SEMIMAJOR_RADIUS, rad2deg(inclination_rad)
                )
            else
                eccentricity, argument_of_perigee = (10^-7, 90)
            end
            push!(
                TLEs,
                TLE(;
                    name=name,
                    epoch_year=25,
                    epoch_day=1,
                    inclination=inclination_rad/π*180,
                    raan=Ω,
                    eccentricity=eccentricity,
                    argument_of_perigee=argument_of_perigee,
                    mean_anomaly=anomaly,
                    mean_motion=SECONDS_PER_DAY*√(
                        G*EARTH_MASS_KG/(altitude_km*1000+SEMIMAJOR_RADIUS)^3
                    )/2/π,
                ),
            )
        end
    end
    return TLEs
end

"""
    generate_regular_array(; <keyword arguments>)

Generate synthetic propagators for a regular constellation of satellites.

# Keyword Arguments

  - `name_prefix::String="SAT"`: Prefix for satellite names.
  - `orbital_planes::Int=12`: Number of orbital planes.
  - `sats_per_plane::Int=18`: Number of satellites per orbital plane.
  - `altitude_km::Int=500`: Altitude of the orbits in kilometers.
  - `inclination_rad::Float64=π/2`: Inclination angle in radians.
  - `frozen_orbits::Bool=false`: Whether to use frozen orbit parameters (see
    [SatelliteAnalysis.jl/frozen_orbits](https://juliaspace.github.io/SatelliteAnalysis.jl/stable/man/frozen_orbits/)).

# Returns

  - `Vector{Propagator}`: Synthetic Propagator objects.
"""
function generate_regular_array(;
    name_prefix::String="SAT",
    orbital_planes::Int=12,
    sats_per_plane::Int=18,
    altitude_km::Int=500,
    inclination_rad::Float64=π/2,
    frozen_orbits::Bool=false,
)
    TLEs = generate_regular_array_TLEs(;
        name_prefix=name_prefix,
        orbital_planes=orbital_planes,
        sats_per_plane=sats_per_plane,
        altitude_km=altitude_km,
        inclination_rad=inclination_rad,
        frozen_orbits=frozen_orbits,
    )
    return [Propagators.init(Val(:SGP4), tle) for tle in TLEs]
end