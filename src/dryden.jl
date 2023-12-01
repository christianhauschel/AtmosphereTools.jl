using Random
using ControlSystems
using Statistics


"""
    dryden(t, h_SI, V_SI;
        W_20_SI=15,
        positive_q::Bool=true,
        positive_r::Bool=true,
        b_SI=1.0,
        seed_values=[23341, 23342, 23343, 23344],
    )

Simulate the Dryden gust model according to MIL-F-8785C [1].

# Arguments
- `t::Vector`: Time vector
- `h::Number`: Altitude [m]
- `V::Number`: Airspeed [m/s]
- `W6`: Wind velocity at 20 ft (≈ 6 m) [m/s] 7.7 light, 15.4 moderate, 23.2 severe
- `positive_q`: Set to `true` if positive q is desired
- `positive_r`: Set to `true` if positive r is desired
- `b_SI`: Wing span [m]
- `seeds::Vector{Int}`: Seed values for the random number generator

# Returns
- `uvw::Matrix`: Linear velocities [m/s]
- `pqr::Matrix`: Angular rates [rad/s]

# References
[1] https://ch.mathworks.com/help/aeroblks/drydenwindturbulencemodelcontinuous.html
"""
function dryden(
    t::Vector,
    h::Number,
    V::Number;
    W6=15.0,
    positive_q::Bool=true,
    positive_r::Bool=true,
    b_SI=1.0,
    seeds::Vector{Int}=[23341, 23342, 23343, 23344],
)
    n = length(t)
    Δt = t[2] - t[1]
    white_noise_u = _whiteNoise(n, Δt, noise_power=π, seed=seeds[1])
    white_noise_v = _whiteNoise(n, Δt, noise_power=π, seed=seeds[2])
    white_noise_w = _whiteNoise(n, Δt, noise_power=π, seed=seeds[3])
    white_noise_r = _whiteNoise(n, Δt, noise_power=π, seed=seeds[4])

    return _dryden(
        t, h, V,
        white_noise_u, white_noise_v, white_noise_w, white_noise_r;
        W_20_SI=W6,
        positive_q=positive_q,
        positive_r=positive_r,
        b_SI=b_SI,
    )

end

function _whiteNoise(n, Δt; noise_power=π, seed=23341)
    Random.seed!(seed)
    noise_multiplier = √noise_power / √Δt
    return randn(n) * noise_multiplier
end


function _dryden(
    t,
    h_SI,
    V_SI,
    white_noise_u::Vector,
    white_noise_v::Vector,
    white_noise_w::Vector,
    white_noise_r::Vector;
    W_20_SI=15,
    positive_q::Bool=true,
    positive_r::Bool=true,
    b_SI=1.0,
)

    white_noise_u = reshape(white_noise_u, 1, length(white_noise_u))
    white_noise_v = reshape(white_noise_v, 1, length(white_noise_v))
    white_noise_w = reshape(white_noise_w, 1, length(white_noise_w))
    white_noise_r = reshape(white_noise_r, 1, length(white_noise_r))


    b = m2ft(b_SI)
    V = ms2fts.(V_SI)
    h = m2ft.(h_SI)
    W_20 = ms2fts.(W_20_SI)

    if h > 1000.0
        error("Altitude is above the low-altitude limit of 304.8 m (1000 ft)!")
    end

    n = length(t)

    # Low-altitude model (< 1000 ft)
    L_w = h
    L_u = L_v = h ./ (0.177 .+ 0.000823 .* h) .^ 1.2

    σ_w = 0.1 .* W_20
    σ_u = σ_w ./ (0.177 .+ 0.000823 .* h) .^ 0.4
    σ_v = σ_w ./ (0.177 .+ 0.000823 .* h) .^ 0.4

    # -------------------------------------
    # Linear velocities
    # -------------------------------------

    # Longitudinal
    H_u = tf([σ_u * sqrt(2 * L_u / π / V)], [L_u / V, 1,])

    # Lateral
    H_v = tf(σ_v * sqrt(L_v / π / V) .* [sqrt(3) * L_v / V, 1], [(L_v / V)^2, 2 * L_v / V, 1])

    # Vertical
    H_w = tf(σ_w * sqrt(L_w / π / V) .* [sqrt(3) * L_w / V, 1], [(L_w / V)^2, 2 * L_w / V, 1])


    # -------------------------------------
    # Angular rates
    # -------------------------------------

    # Longitudinal
    H_p = tf(σ_w * sqrt(0.8 / V) * [(π / (4b))^(1 / 6)], L_w^(1 / 3) * [4b / (π * V), 1])


    # Lateral
    if positive_r
        multiplier_r = 1.0
    else
        multiplier_r = -1.0
    end
    H_r = tf([multiplier_r / V, 0.0], [3 * b / (π * V), 1.0]) #* H_v


    # Vertical
    if positive_q
        multiplier_q = 1.0
    else
        multiplier_q = -1.0
    end
    H_q = tf([multiplier_q / V, 0.0], [4 * b / (π * V), 1.0]) #* H_w


    # -------------------------------------
    # Simulation
    # -------------------------------------

    # Linear
    u, t_u, _x, _u = lsim(H_u, white_noise_u, t)
    v, t_v, _, _ = lsim(H_v, white_noise_v, t)
    w, t_w, _, _ = lsim(H_w, white_noise_w, t)

    u = vec(u)
    v = vec(v)
    w = vec(w)

    u_SI = fts2ms.(u)
    v_SI = fts2ms.(v)
    w_SI = fts2ms.(w)

    # Angular
    p, _, _, _ = lsim(H_p, white_noise_r, t)
    r, _, _, _ = lsim(H_r, reshape(v, 1, n), t)
    q, _, _, _ = lsim(H_q, reshape(w, 1, n), t)
    p = vec(p)
    r = vec(r)
    q = vec(q)

    return u_SI, v_SI, w_SI, p, q, r
end