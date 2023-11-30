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


fname_V_ref = "/mnt/a/Code/10_drones/drone_simulation/dryden_V.csv"
fname_omega_ref = "/mnt/a/Code/10_drones/drone_simulation/dryden_omega.csv"
fname_noise = "/mnt/a/Code/10_drones/drone_simulation/noise.csv"

fs = 1000  # Sampling frequency in Hz
Δt = 1.0 / fs   # Sampling interval
t1 = 0.0
t2 = 5.0
DCM = nothing # trafo matrix: {I} -> {B}

b_SI = 1.0
V_SI = 10.0
h_SI = 5.0
W_20_SI = 15 # wind velocity at 20 ft (≈ 6 m) [m/s] 7.7 light, 15.4 moderate, 23.2 severe

# ==============================================================================
# Calculations
# ==============================================================================

kn2ms(x) = x / 1.94384
ms2kn(x) = x * 1.94384

ft2m(x) = x / 3.28084
m2ft(x) = 3.28084 * x

ms2fts(x) = 3.28084 * x
fts2ms(x) = x / 3.28084

rms(x) = sqrt(mean(x.^2))

# Conversion 

b = m2ft(b_SI)
V = ms2fts(V_SI)
h = m2ft(h_SI)
W_20 = ms2fts(W_20_SI)



if h > 1000.0
    error("Altitude is above the low-altitude limit of 304.8 m (1000 ft)!")
end

# Time
t = Vector(t1:Δt:t2)
n = length(t)

seed_values = [23341, 23342, 23343, 23344] # MATLAB [23341, 23342, 23343, 23344]

# Input
Random.seed!(seed_values[1])
noise_power = π
noise_multiplier = √noise_power / √Δt
white_noise_u = reshape(randn(n), 1, n) * noise_multiplier
Random.seed!(seed_values[2])
white_noise_v = reshape(randn(n), 1, n) * noise_multiplier
Random.seed!(seed_values[3])
white_noise_w = reshape(randn(n), 1, n) * noise_multiplier
Random.seed!(seed_values[4])
white_noise_roll = reshape(randn(n), 1, n) * noise_multiplier

# df_noise = CSV.read(fname_noise, DataFrame)
# df_noise.t = df_noise.time .- df_noise.time[1]

# white_noise_u = df_noise[!, "White Noise:1(1)"]
# white_noise_v = df_noise[!, "White Noise:1(2)"]
# white_noise_w = df_noise[!, "White Noise:1(3)"]
# white_noise_roll = df_noise[!, "White Noise:1(4)"]
# white_noise_u = reshape(white_noise_u, 1, n)
# white_noise_v = reshape(white_noise_v, 1, n)
# white_noise_w = reshape(white_noise_w, 1, n)
white_noise_roll = reshape(white_noise_roll, 1, n)


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
H_p = tf(σ_w * sqrt(0.8 / V) * [(π / (4b))^(1 / 6)], L_w^(1 / 3) * [1, 4b / (π * V)])

Random.seed!(seed_values[4])
random_vector = rand([-1, 1], n)

# Lateral
H_r = tf([0.0, -1 / V], [1.0, 3 * b / (π * V)]) * H_v


# Vertical
H_q = tf([0.0, +1 / V], [1.0, 4 * b / (π * V)]) * H_w



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


# fig1, ax = pplt.subplots(figsize = (7, 3))
# ax[1].plot(t_u, u_SI, lw = 1, label = "u", color="C0")
# ax[1].plot(t_v, v_SI, lw = 1, label = "v", color="C1")
# ax[1].plot(t_w, w_SI, lw = 1, label = "w", color="C2")
# ax[1].plot(t_ref, u_ref, lw = 1, label = "u ref", color="C0", alpha=0.5)
# ax[1].plot(t_ref, v_ref, lw = 1, label = "v ref", color="C1", alpha=0.5)
# ax[1].plot(t_ref, w_ref, lw = 1, label = "w ref", color="C2", alpha=0.5)
# ax[1].set(xlabel = "t [s]", ylabel = "u, v, w [m/s]", title = "Dryden")
# ax[1].legend(ncols = 2)
# fig1

# function calc_Xss(u, fs)
#     X = fft(u)
#     N = length(u)
#     n_freq = div(N, 2) + 1
#     f_nyquist = fs / 2.0
#     f = LinRange(0, f_nyquist, n_freq)
#     Δf = f[2] - f[1]

#     X_ss = X[1:div(N, 2)+1]
#     X_ss *= 2
#     X_ss[1] = X_ss[1] / 2.0
#     return X_ss, f
# end

# U_ss, f = calc_Xss(u, fs)
# U_ss_ref, f_ref = calc_Xss(u_ref, 1.0)


# fig, ax = pplt.subplots(figsize = (7, 3))
# ax[1].plot(f, abs.(U_ss)/length(u), lw = 1, label = "u", color="C0")
# ax[1].plot(f_ref, abs.(U_ss_ref)/length(u_ref), lw = 1, label = "u ref", color="C0", alpha=0.5)
# ax[1].set(
#     xlabel = "f [Hz]",
#     ylabel = "U [m/s]",
#     xscale = "log",
#     # yscale = "log",
#     # xlim = (0.1, 100),
#     # ylim = (1e-4, 1e-1),

# )
# fig


p, t_p, _, _ = lsim(H_p, white_noise_roll, t)
r, t_r, _, _ = lsim(H_r, white_noise_roll, t)
q, t_q, _, _ = lsim(H_q, white_noise_roll, t)

p = vec(p)
r = vec(r)
q = vec(q)

# r = r .* random_vector
# q = q .* random_vector

fig2, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t_p, p, lw=1, label="p", color="C0")
# ax[1].plot(t_r, r, lw=1, label="r", color="C1)
# ax[1].plot(t_q, q, lw=1, label="q", color="C2)
ax[1].plot(t_ref, p_ref, lw=1, label="p ref", alpha=0.5, color="C0")
# ax[1].plot(t_ref, r_ref, lw=1, label="r ref", alpha=0.5, color="C1)
# ax[1].plot(t_ref, q_ref, lw=1, label="q ref", alpha=0.5, color="C2)
ax[1].set(
    xlabel="t [s]",
    ylabel="p, r, q [rad/s]",
)
ax[1].legend(ncols=1)
fig2