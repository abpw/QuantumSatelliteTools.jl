using ..FreespaceChannels: FreespaceChannel

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
Convert a loss value from decibels to probability.

# Arguments
- `loss_dB::Number`: The loss in decibels.

# Returns
- The loss as a probability (0 - 1).
"""
function dB_to_prob(loss_dB::Number)
    1 - 10^(loss_dB / -10)
end

"""
Convert a loss value from probability to decibels.

# Arguments
- `loss_prob::Number`: The loss as a probability (0 - 1).

# Returns
- The loss in decibels.
"""
function prob_to_dB(loss_prob::Number)
    -10 * log10(1 - loss_prob)
end

"""
Calculate the geometric loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total geometric loss.
"""
function geometric_loss_dB(channel::FreespaceChannel)
    20*log10((channel.transmitter_diameter_m + channel.distance_m * waist_radius(channel)) / channel.receiver_diameter_m)
end

"""
Calculate the distance over which the channel passes through the atmosphere.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.
- `atmosphere_height_m`: The altitude of the boundary of the atmosphere in meters.

# Returns
- The distance that the channel passes through the atmosphere in meters.
"""
function atmosphere_distance(channel::FreespaceChannel, atmosphere_height_m::Float64=18e3)
    if channel.elevation_angle_rad == 0
        atmospheric_distance = 0
    else
        relative_elevation = channel.distance_m * sin(channel.elevation_angle_rad)
        atmospheric_elevation = min(relative_elevation, atmosphere_height_m - channel.min_altitude_m)
        atmospheric_distance = atmospheric_elevation / sin(channel.elevation_angle_rad)
    end
    atmospheric_distance
end

"""
Calculate the atmospheric loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total atmospheric loss.
"""
function atmospheric_loss_dB(channel::FreespaceChannel)
    atmospheric_distance = atmosphere_distance(channel)
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
        throw("Molecular absorption is not defined for other wavelengths")
    end

    return molecular_absorption * (atmospheric_distance / 1000)

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

    # transmittance_loss + weather_loss
end

"""
Calculate the pointing loss for a channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.
- `pointing_jitter::Number`: The pointing jitter in microradians (default 5).

# Returns
- The pointing loss in decibels (dB), where larger values mean more loss.
"""
function pointing_loss_dB(channel::FreespaceChannel, pointing_jitter::Float64=5.0)
    # Pointing transmittance factor (0-1)
    η_pnt = exp(-8 * (pointing_jitter * 1e-6)^2 / waist_radius(channel)^2)
    # Convert transmittance to loss in dB so it can be added to other dB losses.
    return -10 * log10(η_pnt)
end

"""
Calculate the reflector loss for a satellite.

# Returns
- The reflector loss for one satellite in decibels.
"""
function reflector_loss()
    # decibels_to_probability(5.854678746311231)
    # 5% loss from reflector paper
    0.05
end

"""
Calculate the swapping loss for a ground station as a probability.

# Arguments
- `memory_read_write_loss::Number`: Read/write loss of the quantum memory as a probability.

# Returns
- The swapping loss for one ground station as a probability.
"""
function swapping_loss(memory_read_write_loss::Float64=0.8)
    # 0.5 is Bell state measurement
    0.5 * memory_read_write_loss
end

"""
Calculate the swapping loss for a ground station in decibels.

# Arguments
- `memory_read_write_loss::Number`: Read/write loss of the quantum memory as a probability.

# Returns
- The swapping loss for one ground station in decibels.
"""
function swapping_loss_dB(memory_read_write_loss::Float64=0.8)
    return prob_to_dB(swapping_loss(memory_read_write_loss))
end

"""
Calculate total loss in a freespace channel.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total loss in decibels.
"""
function total_loss_dB(channel::FreespaceChannel)
    # L_tot = L_geo + L_atm + L_pnt
    g = geometric_loss_dB(channel)
    a = atmospheric_loss_dB(channel)
    p = pointing_loss_dB(channel)

    return g + a + p
end

"""
Calculate total loss in a freespace channel as a probability.

# Arguments
- `channel::FreespaceChannel`: A freespace channel.

# Returns
- The total loss as a probability (0 - 1).
"""
function total_loss(channel::FreespaceChannel)
    dB_to_prob(total_loss_dB(channel))
end