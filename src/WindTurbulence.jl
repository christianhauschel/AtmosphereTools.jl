module WindTurbulence

include("conversions.jl")
export kn2ms, ms2kn, ft2m, m2ft, ms2fts, fts2ms

include("dryden.jl")
export dryden

end # module WindTurbulence
