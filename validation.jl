"""
Dryden model, MIL-F-8785C, low altitude.
# References
https://ch.mathworks.com/help/aeroblks/drydenwindturbulencemodelcontinuous.html
"""


using Random
using FFTW
using ControlSystems
using PyPlot, PyCall
using DataFrames, CSV
using Statistics
using AcousticAnalysis
pplt = pyimport("proplot")
pplt.close("all")

# ==============================================================================
# Parameters
# ==============================================================================


fname_V_ref = "validation/V.csv"
fname_omega_ref = "validation/omega.csv"
fname_noise = "validation/noise.csv"

positive_q = false
positive_r = true

fs = 1000  # Sampling frequency in Hz
Δt = 1.0 / fs   # Sampling interval
t1 = 0.0
t2 = 5.0

b_SI = 1.0
V_SI = 10.0
h_SI = 5.0
W_20_SI = 15 # wind velocity at 20 ft (≈ 6 m) [m/s] 7.7 light, 15.4 moderate, 23.2 severe


df_noise = CSV.read(fname_noise, DataFrame)
white_noise_u = df_noise.u
white_noise_v = df_noise.v
white_noise_w = df_noise.w
white_noise_r = df_noise.r
white_noise_u = reshape(white_noise_u, 1, n)
white_noise_v = reshape(white_noise_v, 1, n)
white_noise_w = reshape(white_noise_w, 1, n)
white_noise_r = reshape(white_noise_r, 1, n)

fig, ax = pplt.subplots(figsize = (7, 8), ncols=1, nrows=4, sharex=true, sharey=false)
ax[1].plot(t, white_noise_u_OG[1, :], lw = 1, label = "Julia", color="C0")
ax[1].plot(t, white_noise_u[1, :], lw = 1, label = "MATLAB", color="C1")
ax[1].set(xlabel = "t [s]", ylabel = "noise u")
ax[2].plot(t, white_noise_v_OG[1, :], lw = 1, label = "Julia", color="C0")
ax[2].plot(t, white_noise_v[1, :], lw = 1, label = "MATLAB", color="C1")
ax[2].set(xlabel = "t [s]", ylabel = "noise v")
ax[3].plot(t, white_noise_w_OG[1, :], lw = 1, label = "Julia", color="C0")
ax[3].plot(t, white_noise_w[1, :], lw = 1, label = "MATLAB", color="C1")
ax[3].set(xlabel = "t [s]", ylabel = "noise w")
ax[4].plot(t, white_noise_r_OG[1, :], lw = 1, label = "Julia", color="C0")
ax[4].plot(t, white_noise_r[1, :], lw = 1, label = "MATLAB", color="C1")
ax[4].set(xlabel = "t [s]", ylabel = "noise r")
ax[1].legend(ncols = 1)
fig

# Low-altitude model (< 1000 ft)
L_w = h
L_u = L_v = h / (0.177 + 0.000823 * h)^1.2

σ_w = 0.1 * W_20
σ_u = σ_w / (0.177 + 0.000823 * h)^0.4
σ_v = σ_w / (0.177 + 0.000823 * h)^0.4


# -------------------------------------
# Linear velocities
# -------------------------------------

# Longitudinal
H_u = tf([σ_u * sqrt(2 * L_u / π / V)], [L_u / V, 1, ])

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



u, t_u, _x, _u = lsim(H_u, white_noise_u, t)
v, t_v, _, _ = lsim(H_v, white_noise_v, t)
w, t_w, _, _ = lsim(H_w, white_noise_w, t)
u = vec(u)
v = vec(v)
w = vec(w)

u_SI = fts2ms.(u)
v_SI = fts2ms.(v)
w_SI = fts2ms.(w)


df_V = CSV.read(fname_V_ref, DataFrame)
u_ref = df_V.u
v_ref = df_V.v
w_ref = df_V.w
df_V.t = convert.(Float64, df_V.t)
t_ref = df_V.t

df_omega = CSV.read(fname_omega_ref, DataFrame)
p_ref = df_omega.p
q_ref = df_omega.q
r_ref = df_omega.r


fig1, ax = pplt.subplots(figsize = (7, 3))
ax[1].plot(t_u, u_SI, lw = 1, label = "u", color="C0")
ax[1].plot(t_v, v_SI, lw = 1, label = "v", color="C1")
ax[1].plot(t_w, w_SI, lw = 1, label = "w", color="C2")
ax[1].plot(t_ref, u_ref, lw = 1, label = "u ref", color="C0", alpha=0.5)
ax[1].plot(t_ref, v_ref, lw = 1, label = "v ref", color="C1", alpha=0.5)
ax[1].plot(t_ref, w_ref, lw = 1, label = "w ref", color="C2", alpha=0.5)
ax[1].set(xlabel = "t [s]", ylabel = "u, v, w [m/s]", title = "Dryden")
ax[1].legend(ncols = 2)
fig1


p, t_p, _, _ = lsim(H_p, white_noise_r, t)
r, t_r, _, _ = lsim(H_r, reshape(v, 1, n), t)
q, t_q, _, _ = lsim(H_q, reshape(w, 1, n), t)

p = vec(p)
r = vec(r)
q = vec(q)

# r = r .* random_vector
# q = q .* random_vector

fig2, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t_p, p, lw=1, label="p", color="C0")
ax[1].plot(t_r, r, lw=1, label="r", color="C1")
ax[1].plot(t_q, q, lw=1, label="q", color="C2")
ax[1].plot(t_ref, p_ref, lw=1, label="p ref", color="C0", ls="--", alpha=0.5)
ax[1].plot(t_ref, r_ref, lw=1, label="r ref", color="C1", ls="--", alpha=0.5)
ax[1].plot(t_ref, q_ref, lw=1, label="q ref", color="C2", ls="--", alpha=0.5)
ax[1].set(
    xlabel="t [s]",
    ylabel="p, r, q [rad/s]",
)
ax[1].legend(ncols=1)
fig2