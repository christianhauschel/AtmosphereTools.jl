"""
    horizontal_wind_model(velocity_abs::Vector, direction::Vector)

Return the wind speed in {I}-frame given the direction and absolute velocity.

# Arguments 
- `velocity_abs::Vector`: absolute velocity of the wind
- `direction::Vector`: direction of the wind [rad], measured from north (y-axis), 
    positive clockwise, direction from which the wind is coming

    Example: Wind coming from north (y-axis) has a direction of 0° and blows from N to S.

# Returns 
- `wind::Matrix`: wind speed in {I} frame, size (n_τ, 3)

# References 
1. https://en.wikipedia.org/wiki/Wind_direction
2. https://ch.mathworks.com/help/aeroblks/horizontalwindmodel.html
"""
function horizontal_wind_model(velocity_abs::Vector, direction::Vector)
    u = -velocity_abs .* sin.(direction)
    v = -velocity_abs .* cos.(direction)
    w = zeros(length(velocity_abs))
    
    # return matrix of size (n_τ, 3)
    return [u v w]
end


"""
    shear(altitude, V_wind_at6m, direction_wind_at6m)

Calculate the wind shear from conditions measured at a height of 6 m (20 ft) above the ground.

# Arguments 
- `altitude`: altitude of the aircraft [m]
- `V_wind_at6m`: wind speed at 6 m (20 ft) [m/s]
- `direction_wind_at6m`: wind direction at 6 m (20 ft) [rad], measured from north (y-axis), 
    positive clockwise, direction from which the wind is coming

    Example: Wind coming from north (y-axis) has a direction of 0° and blows from N to S.

# Keyword Arguments
- `flight_phase::String`: flight phase, "C" (terminal flight phase) or "other" (default: "C")

    Category C flight phases are defined in reference [2] to be terminal flight 
    phases, which include takeoff, approach, and landing.

# Returns
- wind speed at given altitude

# References 
1. https://ch.mathworks.com/help/aeroblks/windshearmodel.html
2. U.S. Military Specification MIL-F-8785C, November 5, 1980.
"""
function shear(altitude, V_wind_at6m; flight_phase = "C")

    W20 = ms2fts(V_wind_at6m)
    h = m2ft(altitude)

    if any(h > 1000.0)
        error("Altitude is above the low-altitude limit of 304.8 m (1000 ft)!")
    end
    if any(h < 3)
        println("WARNING: Altitude is below the low-altitude limit of 3 m (10 ft)!")
    end

    if flight_phase == "C"
        z0 = 0.15
    else
        z0 = 2.0
    end

    u_w = W20 * log(h / z0) / log(20/z0)

    return fts2ms(u_w)
end