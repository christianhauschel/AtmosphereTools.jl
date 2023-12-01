using PyPlot, PyCall
using DataFrames, CSV
using WindTurbulence: dryden
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

u, v, w, p, q, r = dryden(
    t,
    h,
    V;
    W6=W_20_SI,
    positive_q=positive_q,
    positive_r=positive_r,
    b_SI=b_SI,
    seeds=[23341, 23342, 23343, 23344],
)

# fig, ax = pplt.subplots(figsize = (7, 8), ncols=1, nrows=4, sharex=true, sharey=false)
# ax[1].plot(t, white_noise_u_OG[1, :], lw = 1, label = "Julia", color="C0")
# ax[1].plot(t, white_noise_u[1, :], lw = 1, label = "MATLAB", color="C1")
# ax[1].set(xlabel = "t [s]", ylabel = "noise u")
# ax[2].plot(t, white_noise_v_OG[1, :], lw = 1, label = "Julia", color="C0")
# ax[2].plot(t, white_noise_v[1, :], lw = 1, label = "MATLAB", color="C1")
# ax[2].set(xlabel = "t [s]", ylabel = "noise v")
# ax[3].plot(t, white_noise_w_OG[1, :], lw = 1, label = "Julia", color="C0")
# ax[3].plot(t, white_noise_w[1, :], lw = 1, label = "MATLAB", color="C1")
# ax[3].set(xlabel = "t [s]", ylabel = "noise w")
# ax[4].plot(t, white_noise_r_OG[1, :], lw = 1, label = "Julia", color="C0")
# ax[4].plot(t, white_noise_r[1, :], lw = 1, label = "MATLAB", color="C1")
# ax[4].set(xlabel = "t [s]", ylabel = "noise r")
# ax[1].legend(ncols = 1)
# fig




fig1, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t_u, u, lw=1, label="u", color="C0")
ax[1].plot(t_v, v, lw=1, label="v", color="C1")
ax[1].plot(t_w, w, lw=1, label="w", color="C2")
# ax[1].plot(t_ref, u_ref, lw = 1, label = "u ref", color="C0", alpha=0.5)
# ax[1].plot(t_ref, v_ref, lw = 1, label = "v ref", color="C1", alpha=0.5)
# ax[1].plot(t_ref, w_ref, lw = 1, label = "w ref", color="C2", alpha=0.5)
ax[1].set(xlabel="t [s]", ylabel="u, v, w [m/s]", title="Dryden")
ax[1].legend(ncols=2)
fig1



fig2, ax = pplt.subplots(figsize=(7, 3))
ax[1].plot(t_p, p, lw=1, label="p", color="C0")
ax[1].plot(t_r, r, lw=1, label="r", color="C1")
ax[1].plot(t_q, q, lw=1, label="q", color="C2")
# ax[1].plot(t_ref, p_ref, lw=1, label="p ref", color="C0", ls="--", alpha=0.5)
# ax[1].plot(t_ref, r_ref, lw=1, label="r ref", color="C1", ls="--", alpha=0.5)
# ax[1].plot(t_ref, q_ref, lw=1, label="q ref", color="C2", ls="--", alpha=0.5)
ax[1].set(
    xlabel="t [s]",
    ylabel="p, r, q [rad/s]",
)
ax[1].legend(ncols=1)
fig2

# fig1.savefig("validation/dryden_uvw.png", dpi=300)
# fig2.savefig("validation/dryden_prq.png", dpi=300)