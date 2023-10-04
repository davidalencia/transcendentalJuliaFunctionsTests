module ElementalFunctionsTests

export γsectionCC, SSR, DSR
include("random.jl")

export mpfrfunctions, mpfrfunctions_symbol
include("constants.jl")

export ulpdistance, distance2inf, ulperror, ulpdistance32, distance2inf32, ulperror32, mpfrfun
include("zimmermannInterface.jl")

export openinfinite, domain, testRandomBatch, TestsResults,
  save_testRandomBatch, Floating64Distribution, FloatingDistribution, γSection64, γSection
include("testing.jl")





end # module ElementalFunctionsTests
