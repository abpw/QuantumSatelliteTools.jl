using GeoDatasets: landseamask

############################################################################################
#                                     Earth Constants                                      #
############################################################################################

"""
    AU_M

Astronomical unit [m].
"""
const AU_M = 1.495978707e11

"""
    SEMIMAJOR_RADIUS

Semi-major axis of Earth [m].
"""
const SEMIMAJOR_RADIUS = 6378137.0

"""
    SEMIMINOR_RADIUS

Semi-minor axis of Earth [m].
"""
const SEMIMINOR_RADIUS = 6356752.3142

"""
    EARTH_MASS_KG

Mass of Earth [kg].
"""
const EARTH_MASS_KG = 5.9722e24

"""
    EQUATORIAL_CIRCUMFERENCE_KM

Equatorial circumference of Earth [km].
"""
const EQUATORIAL_CIRCUMFERENCE_KM = 2π * SEMIMAJOR_RADIUS / 1000

############################################################################################
#                                     Other Constants                                      #
############################################################################################

"""
    G

Gravitational constant [m³ kg⁻¹ s⁻²].
"""
const G = 6.6743e-11

"""
    SECONDS_PER_DAY

Number of seconds in a day.
"""
const SECONDS_PER_DAY = 60 * 60 * 24

"""
    SIN_60

Sine of 60 degrees.
"""
const SIN_60 = √3 / 2

"""
    LAND_SEA_MASK

Land-sea mask dataset from [GeoDatasets.jl](https://github.com/JuliaGeo/GeoDatasets.jl),
used to determine if a geographic coordinate falls on land or water. The mask is a 2D array
with structure `[[lats], [lons], [data]]`. A `data` point is `0` for ocean, `1` for land, or
`2` for lake.
"""
const LAND_SEA_MASK = landseamask(; resolution=('f'), grid=1.25)