# PASSED TESTS

function dist(city_1::NTuple{2,Float64}, city_2::NTuple{2,Float64})
    return have_R_sine_distance(city_1[1],city_1[2], city_2[1], city_2[2])
end


function have_R_sine_distance(lat1, lon1, lat2, lon2; R::Float64=6_371_000.0)
    lat1r = deg2rad(Float64(lat1))
    lon1r = deg2rad(Float64(lon1))
    lat2r = deg2rad(Float64(lat2))
    lon2r = deg2rad(Float64(lon2))

    dlat = lat2r - lat1r
    dlon = lon2r - lon1r

    s = sin(dlat/2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon/2)^2
    c = 2 * asin(min(1.0, sqrt(s)))  # clamp for numerical safety
    return R * c
end