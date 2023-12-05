using PyPlot, PyCall
using DataFrames, CSV
using WindTurbulence: dryden, _dryden
pplt = pyimport("proplot")
pplt.close("all")

# ==============================================================================
# Parameters
# ==============================================================================

fname_V_ref = "validation/V.csv"
fname_omega_ref = "validation/omega.csv"
fname_noise = "validation/noise.csv"


fs = 1000  # Sampling frequency in Hz
Δt = 1.0 / fs   # Sampling interval
t1 = 0.0
t2 = 5.0
t = Vector(t1:Δt:t2)
n = length(t)

unit = ones(n)

b = 1.0
V = 10.0unit
h = 5.0unit
velocity_6m = 15 # wind velocity at 20 ft (≈ 6 m) [m/s] 7.7 light, 15.4 moderate, 23.2 severe


df_noise = CSV.read(fname_noise, DataFrame)
noise_power = π
noise_multiplier_MATLAB = √noise_power / √Δt
noise_u_ref = df_noise.u / noise_multiplier_MATLAB
noise_v_ref = df_noise.v / noise_multiplier_MATLAB
noise_w_ref = df_noise.w / noise_multiplier_MATLAB
noise_p_ref = df_noise.r / noise_multiplier_MATLAB

df_Vref = CSV.read(fname_V_ref, DataFrame)
u_ref = df_Vref.u
v_ref = df_Vref.v
w_ref = df_Vref.w

df_omegaref = CSV.read(fname_omega_ref, DataFrame)
p_ref = df_omegaref.p
q_ref = df_omegaref.q
r_ref = df_omegaref.r


u, v, w, p, q, r = _dryden(
    t, h, V,
    noise_u_ref, noise_v_ref, noise_w_ref, noise_p_ref;
    velocity_6m=velocity_6m,
    wingspan=b,
)

u_test, v_test, w_test, p_test, q_test, r_test = dryden(
    t, h, V,;
    velocity_6m=velocity_6m,
    wingspan=b,
)


fig1, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t, u, lw=1, label="u", color="C0")
ax[1].plot(t, v, lw=1, label="v", color="C1")
ax[1].plot(t, w, lw=1, label="w", color="C2")
ax[1].plot(t, u_ref, lw = 1, label = "u ref", color="C0", alpha=0.5)
ax[1].plot(t, v_ref, lw = 1, label = "v ref", color="C1", alpha=0.5)
ax[1].plot(t, w_ref, lw = 1, label = "w ref", color="C2", alpha=0.5)
ax[1].set(xlabel="t [s]", ylabel="u, v, w [m/s]", title="Dryden")
ax[1].legend(ncols=3)
fig1


fig2, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t, p, lw=1, label="p", color="C0")
ax[1].plot(t, r, lw=1, label="r", color="C1")
ax[1].plot(t, q, lw=1, label="q", color="C2")
ax[1].plot(t, p_ref, lw=1, label="p ref", color="C0", ls="--", alpha=0.5)
ax[1].plot(t, r_ref, lw=1, label="r ref", color="C1", ls="--", alpha=0.5)
ax[1].plot(t, q_ref, lw=1, label="q ref", color="C2", ls="--", alpha=0.5)
ax[1].set(
    xlabel="t [s]",
    ylabel="p, r, q [rad/s]",
)
ax[1].legend(ncols=3)
fig2

fig1.savefig("validation/dryden_uvw.png", dpi=300)
fig2.savefig("validation/dryden_prq.png", dpi=300)


fig3, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t, u_test, lw=1, label="u", color="C0")
ax[1].plot(t, v_test, lw=1, label="v", color="C1")
ax[1].plot(t, w_test, lw=1, label="w", color="C2")
ax[1].set(xlabel="t [s]", ylabel="u, v, w [m/s]", title="Dryden")
ax[1].legend(ncols=3)
fig3

fig4, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t, p_test, lw=1, label="p", color="C0")
ax[1].plot(t, r_test, lw=1, label="r", color="C1")
ax[1].plot(t, q_test, lw=1, label="q", color="C2")
ax[1].set(
    xlabel="t [s]",
    ylabel="p, r, q [rad/s]",
)
ax[1].legend(ncols=3)
fig4