module ElementalFunctionsTests

export @ulperror

export γsectionCC, SSR, DSR
include("random.jl")

export mpfrfunctions
include("constants.jl")

export ulpdistance, distance2inf, ulperror, mpfrfun
include("zimmermannInterface.jl")

export openinfinite, domain, testRandomBatch, save_testRandomBatch, testRandomBatch_DSR, testRandomBatch_γ
include("testing.jl")


# for f in mpfrfunctions[14:length(mpfrfunctions)]
#   save_testRandomBatch(Symbol(f), rand=() -> DSR(domain(Symbol(f))...), imagesDir="results/DSR/")
# end

# testRandomBatch_(Symbol(f), rand=() -> DSR(domain(Symbol(f))...), imagesDir=".", batchSize=1_000_000)
# @time testRandomBatch_γ(cos, openinfinite, batchSize=1)

end # module ElementalFunctionsTests
