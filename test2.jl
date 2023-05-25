using ProgressMeter

const MAX_INT = 2_139_095_040 - 1
# const COMPILER = "/opt/homebrew/opt/llvm@16/bin/clang"
# # es posible encontrar las bibliotecas usando `pkg-config  --cflags --libs mpfr`
# const FLAGS = [
#     "-Wall",
#     "-shared",
#     "-O3",
#     "-o", "check.so",
#     "-lmpfr",
#     "-lgmp",
#     "-lm",
#     "-fopenmp",
#     "-I/opt/homebrew/Cellar/gmp/6.2.1_1/include",
#     "-L/opt/homebrew/Cellar/gmp/6.2.1_1/lib",
#     "-I/opt/homebrew/Cellar/mpfr/4.2.0-p4/include",
#     "-L/opt/homebrew/Cellar/mpfr/4.2.0-p4/lib",
# ]
# const FILE = "check_exhaustive_wrapper_as_argument.c"
const limits = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 8, 16, 32, 64, Inf]
const buckets = zeros(Int, 21)
const FN = ARGS[1]
const fndomain = length(ARGS) >= 2 ? (x) -> x ∈ eval(Meta.parse(ARGS[2])) : (x) -> true



const c_fun = @cfunction(getfield(Main, Symbol(FN)), Cdouble, (Cdouble,))

function check(x::UInt32, fn::String)::Float64
    return ccall(
        (:check, "check2.so"),
        Cdouble,
        (Cuint, Ptr{Cdouble},),
        x, c_fun
    )
end

function saveInBucket(err::Float64)
    if err <= 0.5
        return
    end
    ix = 1
    while err > limits[ix]
        ix += 1
    end
    buckets[ix] += 1
end

function cleanBuckets()
    for i in 1:18
        buckets[i] = 0
    end
end

function saveerrors(fun; domain=(x) -> true)
    maxErr = (-1, 0)
    @showprogress for i in 0:(MAX_INT-1)
        x = convert(UInt32, i)
        if domain(reinterpret(Float32, x))
            err = check(x, fun)
            if (err > maxErr[2])
                maxErr = (i, err)
            end
            saveInBucket(err)
        end
        #negative values
        x = x + 0x80000000
        if domain(reinterpret(Float32, x))
            err = check(x, fun)
            if (err > maxErr[2])
                maxErr = (i, err)
            end
            saveInBucket(err)
        end
    end
    println(fun)
    println("max error: $maxErr")
    println(limits)
    println(buckets)
end

# saveerrors(
#     FN,
#     domain=fndomain
# )
check(convert(UInt32, 5000000), FN) # FN=csc