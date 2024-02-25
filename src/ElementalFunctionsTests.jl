module ElementalFunctionsTests

export γsectionCC, SSR, DSR
include("random.jl")

export mpfrfunctions, mpfrfunctions_symbol, mpfrfunctions_str
include("constants.jl")

export @ulperror, ulpdistance, distance2inf, ulperror, ulpdistance32, distance2inf32, ulperror32, mpfrfun
include("zimmermannInterface.jl")

export Error, Error32, buckets, add2bucket!, openinfinite, domain, testRandomBatch, TestsResults,
  save_testRandomBatch, Floating64Distribution, FloatingDistribution, γSection64, γSection, mean, variance
include("testing.jl")





end # module ElementalFunctionsTests
