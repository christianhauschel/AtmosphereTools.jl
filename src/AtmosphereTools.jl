module AtmosphereTools

include("conversions.jl")
export kn2ms, ms2kn, ft2m, m2ft, ms2fts, fts2ms

include("dryden.jl")
export dryden

include("wind.jl")
export horizontal_wind_model, shear

end # module AtmosphereTools
