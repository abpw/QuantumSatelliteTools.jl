"""
Convert decibel loss to probability (linear) loss.

# Arguments
- `loss_db::Number`: Loss in decibels.

# Returns
- Loss as a probability (0 to 1).
"""
function decibel_to_probability(loss_db::Number)
    1 - 10^(-loss_db / 10)
end

"""
Calculate the beam waist radius for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The half beam divergence (in radians).
"""
function waist_radius(channel::FreespaceChannel)
    # w₀ =  122λ/Dₜ
    1.22 * (channel.wavelength_nm * 1e-9) / channel.transmitter_diameter_m
end

"""
Calculate the geometric loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total geometric loss in decibels.
"""
function geometric_loss(channel::FreespaceChannel)
    20 * log10((channel.transmitter_diameter_m + channel.distance_m * waist_radius(channel)) / channel.receiver_diameter_m)
end

"""
Calculate the distance over which the channel passes through the atmosphere.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.
- `atmosphere_height_m`: The altitude of the boundary of the atmosphere in meters.

# Returns
- The distance that the channel passes through the atmosphere in meters.
"""
function atmosphere_distance(channel::FreespaceChannel, atmosphere_height_m::Num64=18e3)
    if channel.elevation_angle_rad == 0
        # Horizontal link - no meaningful atmosphere distance calculation
        return 0.0
    end
    
    # If satellite is already below atmosphere boundary, signal travels through entire atmosphere
    if channel.min_altitude_m <= atmosphere_height_m
        # Calculate distance from satellite altitude to ground through atmosphere
        altitude_through_atmosphere = atmosphere_height_m
        return altitude_through_atmosphere / sin(channel.elevation_angle_rad)
    else
        # Satellite is above atmosphere - calculate distance through atmosphere only
        # The signal enters atmosphere at atmosphere_height_m and travels to ground
        return atmosphere_height_m / sin(channel.elevation_angle_rad)
    end
end

"""
Calculate the atmospheric loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total atmospheric loss in decibels.
"""
function atmospheric_loss(channel::FreespaceChannel)
    atm_distance = atmosphere_distance(channel)
    # Transmittance loss
    if channel.wavelength_nm == 550
        molecular_absorption = 0.13
    elseif channel.wavelength_nm == 690
        molecular_absorption = 0.01
    elseif channel.wavelength_nm == 850
        molecular_absorption = 0.41
    elseif channel.wavelength_nm == 1550
        molecular_absorption = 0.01
    else
        throw(ArgumentError("Molecular absorption is not defined for wavelength $(channel.wavelength_nm)nm. Supported wavelengths: 550, 690, 850, 1550 nm"))
    end

    transmittance_loss = molecular_absorption * (atm_distance / 1000)

    # Weather loss - to be finished later
    # if channel.conditions == fog
    #     if atmosphere_distance > 50
    #         p = 1.6
    #     elseif atmosphere_distance > 6 && atmosphere_distance < 50
    #         p = 1.3
    #     elseif atmosphere_distance < 6
    #         p = 0.585 * atmosphere_distance^(1 / 3)
    #     else
    #         throw("Fog loss is not defined for distances 6 km or 50 km")
    #     end
    #     weather_loss = (3.91 / atmosphere_distance)(channel.wavelength_nm / 550)^-p
    #     # TODO finish weather
    # elseif channel.conditions == rain
    # elseif channel.conditions == snow
    # else
    #     weather_loss = 0
    # end

    transmittance_loss
    # TODO: Add weather loss when implemented
    # transmittance_loss + weather_loss
end

"""
Calculate the pointing loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.
- `pointing_jitter::Num64`: The pointing jitter in microradians (default 5).

# Returns
- The total pointing loss in decibels.
"""
function pointing_loss(channel::FreespaceChannel, pointing_jitter::Num64=5)
    # L_PNT = exp(-8Θ²ⱼ/w²₀)
    exp(-8 * (pointing_jitter * 10^-6)^2 / waist_radius(channel)^2)

end

"""
Calculate the reflector loss for a satellite.

# Returns
- The reflector loss for one satellite in decibels.
"""
function reflector_loss()
    5.854678746311231
end

"""
Calculate the swapping loss for a ground station.

# Arguments
- `memory_read_write_loss::Num64`: Read/write loss of the quantum memory as a probability.

# Returns
- The swapping loss for one ground station as a probability.
"""
function swapping_loss(memory_read_write_loss::Num64=0.8)
    # 0.5 is Bell state measurement
    0.5 * memory_read_write_loss
end
