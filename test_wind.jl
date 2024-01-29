using Revise
using AtmosphereTools 
using PyPlot, PyCall
pplt = pyimport("proplot")
pplt.close("all")

shear(300, 10)

h = 0:0.1:30
fig, ax = pplt.subplots(figsize=(6,5))
ax[1].plot(h, shear.(h, 7), label="7 m/s")
ax[1].plot(h, shear.(h, 15), label="15 m/s")
ax[1].plot(h, shear.(h, 30), label="30 m/s")
ax[1].legend()
ax[1].format(
    xlabel="Altitude [m]",
    ylabel="Wind speed [m/s]",
    title="Wind shear (MIL-F-8785C)",
)
fig 
# fig.savefig("doc/img/shear.png", dpi=300)