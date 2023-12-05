using Random


"""
    dryden(t::Vector, h::Vector, V::Vector;
        velocity_6m=15.0,
        wingspan=1.0,
        seeds::Vector{Int}=[23341, 23342, 23343, 23344],
    )

Generate a dryden turbulence model for the given time vector `t`, altitude 
vector `h`, and velocity vector `V`.

# Arguments
- `t::Vector`: Time vector (s)
- `h::Vector`: Altitude vector (m)
- `V::Vector`: Aircraft velocity vector (m/s)

# Optional Arguments
- `velocity_6m::Float64`: Velocity at 6 m altitude (m/s), this sets the turbulenc 
    intensity (light turbulence: 7.72 m/s, moderate: 15.4 m/s, severe: 23.2 m/s)
- `wingspan::Float64`: Aircraft wingspan (m)
- `seeds::Vector{Int}`: Random number generator seeds for input noise
"""
function dryden(
    t::Vector,
    h::Vector,
    V::Vector;
    velocity_6m=15.4,
    wingspan=1.0,
    seeds::Vector{Int}=[23341, 23342, 23343, 23344],
)
    n = length(t)

    # White noise input
    Random.seed!(seeds[1])
    noise_u = randn(n)
    Random.seed!(seeds[2])
    noise_v = randn(n)
    Random.seed!(seeds[3])
    noise_w = randn(n)
    Random.seed!(seeds[4])
    noise_r = randn(n)

    return _dryden(
        t, h, V,
        noise_u, noise_v, noise_w, noise_r;
        velocity_6m=velocity_6m,
        wingspan=wingspan,
    )
end

"""
    _dryden(t::Vector, h::Vector, V::Vector
        noise_u::Vector, noise_v::Vector, noise_w::Vector, noise_p::Vector;
        velocity_6m=15.0,
        wingspan=1.0,
        seeds::Vector{Int}=[23341, 23342, 23343, 23344],
    )

Generate a dryden turbulence model for the given time vector `t`, altitude 
vector `h`, and velocity vector `V`.

# Arguments
- `t::Vector`: Time vector (s)
- `h::Vector`: Altitude vector (m)
- `V::Vector`: Aircraft velocity vector (m/s)
- `noise_u::Vector`: White noise input for u (m/s)
- `noise_v::Vector`: White noise input for v (m/s)
- `noise_w::Vector`: White noise input for w (m/s)
- `noise_p::Vector`: White noise input for p (rad/s)

# Optional Arguments
- `velocity_6m::Float64`: Velocity at 6 m altitude (m/s), this sets the turbulenc 
    intensity (light turbulence: 7.72 m/s, moderate: 15.4 m/s, severe: 23.2 m/s)
- `wingspan::Float64`: Aircraft wingspan (m)
- `seeds::Vector{Int}`: Random number generator seeds for input noise

# Returns 
- `u::Vector`: u velocity (m/s)
- `v::Vector`: v velocity (m/s)
- `w::Vector`: w velocity (m/s)
- `p::Vector`: p angular velocity (rad/s)
- `q::Vector`: q angular velocity (rad/s)
- `r::Vector`: r angular velocity (rad/s)
"""
function _dryden(
    t::Vector,
    h_SI::Vector,
    V_SI::Vector,
    noise_u::Vector,
    noise_v::Vector,
    noise_w::Vector,
    noise_p::Vector;
    velocity_6m=15,
    wingspan=1.0,
)
    n = length(t)
    unit = ones(n)

    # Convert to imperial units
    b = m2ft(wingspan)
    V = ms2fts.(V_SI)
    h = m2ft.(h_SI)
    W_20 = ms2fts.(velocity_6m)

    # Check if any element of h > 1000.0
    if any(h .> 1000.0)
        error("Altitude is above the low-altitude limit of 304.8 m (1000 ft)!")
    end

    # Low-altitude model (< 1000 ft)
    L_w = h
    L_u = L_v = h ./ (0.177 .+ 0.000823 .* h) .^ 1.2
    L_p = sqrt.(L_w .* b) ./ 2.6
    σ_w = 0.1 .* W_20 .* unit
    σ_u = σ_w ./ (0.177 .+ 0.000823 .* h) .^ 0.4
    σ_v = σ_w ./ (0.177 .+ 0.000823 .* h) .^ 0.4


    u = zeros(n)
    v = zeros(n)
    w = zeros(n)
    r = zeros(n)
    q = zeros(n)
    p = zeros(n)

    τ_u = L_u ./ V
    τ_v = L_v ./ V
    τ_w = L_w ./ V
    τ_p = L_p ./ V
    τ_q = 4.0b ./ (π .* V)
    τ_r = 3.0b ./ (π .* V)

    σ_r = sqrt.(2.0π ./ (3.0 .* L_w .* b)) .* σ_w
    σ_p = 1.9 ./ sqrt.(L_w .* b) .* σ_w
    σ_q = sqrt.(π ./ (2.0L_w .* b)) .* σ_w

    for i = 2:n
        Δt = t[i] - t[i-1]

        # Linear
        u[i] = (1.0 - Δt / τ_u[i]) * u[i-1] + σ_u[i] * sqrt(2Δt / τ_u[i]) * noise_u[i-1]
        v[i] = (1.0 - Δt / τ_v[i]) * v[i-1] + σ_v[i] * sqrt(2Δt / τ_v[i]) * noise_v[i-1]
        w[i] = (1.0 - Δt / τ_w[i]) * w[i-1] + σ_w[i] * sqrt(2Δt / τ_w[i]) * noise_w[i-1]

        # Rotational
        p[i] = (1.0 - Δt / τ_p[i]) * p[i-1] + σ_p[i] * sqrt(2Δt / τ_p[i]) * noise_p[i]
        r[i] = (1.0 - Δt / τ_r[i]) * r[i-1] + π / (3.0b) * (v[i] - v[i-1])
        q[i] = (1.0 - Δt / τ_q[i]) * q[i-1] + π / (4.0b) * (w[i] - w[i-1])
    end
    u_SI = fts2ms.(u)
    v_SI = fts2ms.(v)
    w_SI = fts2ms.(w)

    return u_SI, v_SI, w_SI, p, q, r
end