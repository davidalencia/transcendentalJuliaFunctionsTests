using Pkg
Pkg.activate("..")
using ElementalFunctionsTests

using ProgressMeter
using JSON
using Printf

JSON.lower(e::Error) =Dict(:x=>@sprintf("%a", e.x), :err=>e.err)

function savearray(arr, file="Float64Errors_GammaSection.json")
    open(file,"w") do f
      write(f, JSON.json(arr))
    end
end

errors = []
for f in mpfrfunctions_symbol
    push!(errors, testRandomBatch(f, γSection(domain(f)...), batchSize=1_000_000_000))
    savearray(errors)
end


