using DifferentialEquations
using PyPlot, PyCall
using DataInterpolations
# using Interpolations
import Gradient
pplt = pyimport("proplot")
pplt.close("all")

tend = 1.0
fs = 50
Δτ = 1.0 / fs
τ = Vector(0.0:Δτ:tend)
n_τ = length(τ)

Random.seed!(12131312)
noise_multiplier = √π / √Δτ
input = randn(n_τ) * noise_multiplier

# H_u
# 14.574876955727152
# ------------------------
# 3.656741260003448s + 1.0
Hu_b0 = 14.57 * ones(n_τ)
Hu_a0 = 1.0 * ones(n_τ)
Hu_a1 = 3.656741260003448 * ones(n_τ)

# H_v
# 65.27468108062853s + 10.305994330354213
# ------------------------------------------------
# 13.371756642611604s^2 + 7.313482520006896s + 1.0
Hv_b0 = 10.305994330354213 * ones(n_τ)
Hv_b1 = 65.27468108062853 * ones(n_τ)
Hv_a0 = 1.0 * ones(n_τ)
Hv_a1 = 7.313482520006896 * ones(n_τ)
Hv_a2 = 13.371756642611604 * ones(n_τ)

# H_r
# 1.9895722156712468s^2 + 0.31412669713714214s
# ------------------------------------------------------------------------
# 1.2769087004961142s^3 + 14.070142779176766s^2 + 7.408975485862033s + 1.0
Hp_b0 = 0.0 * ones(n_τ)
Hp_b1 = 0.31412669713714214 * ones(n_τ)
Hp_b2 = 1.9895722156712468 * ones(n_τ)
Hp_a0 = 1.0 * ones(n_τ)
Hp_a1 = 7.408975485862033 * ones(n_τ)
Hp_a2 = 14.070142779176766 * ones(n_τ)
Hp_a3 = 1.2769087004961142 * ones(n_τ)

function itp(x, y)
    return LinearInterpolation(y, x)
end

function interpolate(x, y)

end

# function solve_LTV(
#     a0::Vector{Float64}, 
#     a1::Vector{Float64}, 
#     b0::Vector{Float64}, 
#     input::Vector{Float64}, 
#     τ::Vector{Float64}; 
#     y0::Float64=0.0
# )
#     f_a0 = itp(τ, a0)
#     f_a1 = itp(τ, a1)
#     f_b0 = itp(τ, b0)
#     f_input = itp(τ, input)

function solve_LTV(
    a::Vector{Vector{Float64}},
    b::Vector{Vector{Float64}},
    input::Vector{Float64},
    τ::Vector{Float64};
    y0::Union{Vector,Nothing}=nothing
)
    m = length(a)
    n = length(b)
    if m <= n
        error("The number of elements in a must be greater than the number of elements in b!")
    end

    vec_fa = Vector(undef, m)
    vec_fb = Vector(undef, n)

    for i = 1:m
        vec_fa[i] = itp(τ, a[i])
    end
    for i = 1:n
        vec_fb[i] = itp(τ, b[i])
    end

    f_input = itp(τ, input)

    # Order of ODE system
    order = m - 1

    if y0 === nothing
        y0 = zeros(order)
    end


    function ode(du, u, p, t)
        if order >= 2
            for i = 1:(order-1)
                du[i] = u[i+1]
            end
        end

        df = input
        for i = 1:n
            du[order] += vec_fb[i](t) * itp(τ, df)(t)
            df = Gradient.gradient(τ, df)
        end
        for i = 1:(m-1)
            du[order] -= vec_fa[i](t) * u[i]
        end
        du[order] /= vec_fa[m](t)

        return du
    end

    tspan = (0.0, tend)
    prob = ODEProblem(ode, y0, tspan)


    sol = solve(prob, saveat=Δτ, adaptive=false, dt=Δτ)
    
    n_τ = length(τ)
    U = zeros(n_τ, order)
    for i = 1:n_τ
        U[i,:] = sol.u[i]
    end
    return U, sol
end

@time Hu_u, Hu_sol = solve_LTV([Hu_a0, Hu_a1], [Hu_b0], input, τ)
@time Hv_u, Hv_sol = solve_LTV([Hv_a0, Hv_a1, Hv_a2], [Hv_b0, Hv_b1], input, τ)
@time Hp_u, Hp_sol = solve_LTV([Hp_a0, Hp_a1, Hp_a2, Hp_a3], [Hp_b0, Hp_b1, Hp_b2], input, τ)


fig, ax = pplt.subplots(figsize=(5, 3))
ax[1].plot(τ, Hu_u[:,1], lw=1)
ax[1].plot(τ, Hv_u[:,1], lw=1)
ax[1].plot(τ, Hp_u[:,1], lw=1)
fig