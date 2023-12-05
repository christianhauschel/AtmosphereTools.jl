using Random

function dryden(
    t::Vector,
    h::Vector,
    V::Vector;
    velocity_6m=15.0,
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