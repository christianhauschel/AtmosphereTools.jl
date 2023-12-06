using Revise
using AtmosphereTools 
using PyPlot, PyCall
pplt = pyimport("proplot")
pplt.close("all")

shear(300, 10)

h = 0:10:300
fig, ax = pplt.subplots(figsize=(6,5))
ax[1].plot(shear.(h, 7), h, label="7 m/s")
ax[1].plot(shear.(h, 15), h, label="15 m/s")
ax[1].plot(shear.(h, 30), h, label="30 m/s")
ax[1].legend()
ax[1].format(
    ylabel="Altitude [m]",
    xlabel="Wind speed [m/s]",
    title="Wind shear (MIL-F-8785C)",
)
fig 
fig.savefig("doc/img/shear.png", dpi=300)