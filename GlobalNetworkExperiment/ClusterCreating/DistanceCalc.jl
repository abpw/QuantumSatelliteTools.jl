function dist(city_1::NTuple{2,Float32}, city_2::NTuple{2,Float64})
    return have_R_sine(city_1[1],city_1[2], city_2[1], city_2[2])
end


function have_R_sine(lat1, lon1, lat2, lon2; R=6371000.0)
    d1, l1, d2, l2 = deg2rad.((lat1, lon1, lat2, lon2))
    ddist, ldist = d2-d1, l2-l1
    a = sin(ddist/2)^2 + cos(d1)*cos(d2)*sin(ldist/2)^2
    return 2R * asin(sqrt(a))
end