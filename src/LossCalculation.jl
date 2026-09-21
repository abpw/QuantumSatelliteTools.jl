using ..FreespaceChannels: FreespaceChannel
using ..helpers: dB_to_prob, prob_to_dB

"""
    waist_radius(channel::FreespaceChannel) -> Float64

Calculate beam waist radius [m] for a `channel`.
"""
function waist_radius(channel::FreespaceChannel)
    # w₀ =  122λ/Dₜ
    return 1.22 * (channel.wavelength_nm * 1e-9) / channel.transmitter_diameter_m
end

"""
    geometric_loss_dB(channel::FreespaceChannel) -> Float64

Calculate geometric loss [dB] for a `channel`.
"""
function geometric_loss_dB(channel::FreespaceChannel)
    return 20*log10(
        (channel.transmitter_diameter_m + channel.distance_m * waist_radius(channel)) /
            channel.receiver_diameter_m,
    )
end

"""
    atmosphere_distance(
        channel::FreespaceChannel;
        atmosphere_height_m::Float64=18e3
    ) -> Float64

Calculate the distance [m] over which the `channel` passes through the atmosphere, which is
bounded by `atmosphere_height_m` [m].
"""
function atmosphere_distance(channel::FreespaceChannel; atmosphere_height_m::Float64=18e3)
    if channel.elevation_angle_rad == 0
        atm_dist = 0
    else
        relative_elev = channel.distance_m * sin(channel.elevation_angle_rad)
        atm_elev = min(
            relative_elev, atmosphere_height_m - channel.min_altitude_m
        )
        atm_dist = atm_elev / sin(channel.elevation_angle_rad)
    end
    return atm_dist
end

"""
    atmospheric_loss_dB(channel::FreespaceChannel) -> Float64

Calculate atmospheric loss [dB] for a `channel`.

# Throws

  - `ErrorException`: If the wavelength is not one of the defined values.
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
        throw(error("Molecular absorption is not defined for other wavelengths"))
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
    pointing_loss_dB(channel::FreespaceChannel; pointing_jitter::Real=2.0) -> Float64

Calculate pointing loss [dB] for a `channel`, with a `pointing_jitter` [μrad].
"""
function pointing_loss_dB(channel::FreespaceChannel; pointing_jitter::Real=2.0)
    # Pointing transmittance factor (0-1)
    η_pnt = exp(-8 * (pointing_jitter * 1e-6)^2 / waist_radius(channel)^2)
    # Convert transmittance to loss in dB so it can be added to other dB losses.
    return -10 * log10(η_pnt)
end

"""
    reflector_loss() -> Float64

Calculate reflector loss [prob] for a satellite.
"""
function reflector_loss()
    # 5% loss from reflector paper
    return 0.05
end

"""
    swapping_loss(; memory_read_write_loss::Real=0.8) -> Float64

Calculate entanglement swapping loss [prob] with a `memory_read_write_loss` [prob].
"""
function swapping_loss(; memory_read_write_loss::Real=0.8)
    # 0.5 is Bell state measurement
    return 0.5 * memory_read_write_loss
end

"""
    swapping_loss_dB(;memory_read_write_loss::Real=0.8) -> Float64

Calculate entanglement swapping loss [dB] with a `memory_read_write_loss` [prob].
"""
function swapping_loss_dB(; memory_read_write_loss::Real=0.8)
    return prob_to_dB(swapping_loss(; memory_read_write_loss=memory_read_write_loss))
end

"""
    total_loss_dB(channel::FreespaceChannel) -> Float64

Calculate total loss [dB] in a freespace `channel`.
"""
function total_loss_dB(channel::FreespaceChannel)
    # L_tot = L_geo + L_atm + L_pnt
    g = geometric_loss_dB(channel)
    a = atmospheric_loss_dB(channel)
    p = pointing_loss_dB(channel)

    return g + a + p
end

"""
    total_loss(channel::FreespaceChannel) -> Float64

Calculate total loss [prob] in a freespace `channel`.
"""
function total_loss(channel::FreespaceChannel)
    return dB_to_prob(total_loss_dB(channel))
end