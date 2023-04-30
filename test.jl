using ProgressMeter
using Base.Threads
using Avro
using DataFrames

const MAX_INT = 2_139_095_040 - 1
const COMPILER = "/opt/homebrew/opt/llvm@16/bin/clang"
# es posible encontrar las bibliotecas usando `pkg-config  --cflags --libs mpfr`
const FLAGS = [
    "-Wall",
    "-shared",
    "-O3",
    "-o", "check.so",
    "-lmpfr",
    "-lgmp",
    "-lm",
    "-fopenmp",
    "-I/opt/homebrew/Cellar/gmp/6.2.1_1/include",
    "-L/opt/homebrew/Cellar/gmp/6.2.1_1/lib",
    "-I/opt/homebrew/Cellar/mpfr/4.2.0-p4/include",
    "-L/opt/homebrew/Cellar/mpfr/4.2.0-p4/lib",
]
const FILE = "check_exhaustive_wrapper_as_argument.c"
const limits = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, Inf]
const buckets = zeros(Int, 17)
# const FN = ARGS[1]


struct UlpError
    reinterpret::Float32
    error::Float64
end

function compile(fn)
    cmd = Cmd([COMPILER, FLAGS..., FILE, "-DSTR=$fn"])
    run(cmd)
end
@generated function createcfun(::Val{fn}) where {fn}
    return :(@cfunction($fn, Cdouble, (Cdouble,)))
end

function check(x::UInt32, fn::String)::Float64
    return ccall(
        (:check, "check.so"),
        Cdouble,
        (Cuint, Ptr{Cdouble},),
        x, createcfun(Val(Symbol(fn)))
    )
end

function getMaxError(fn)
    errmax::Float64 = 0
    err::Float64 = 0
    @showprogress for i in 0:(2_139_095_040-1)
        err = check(convert(UInt32, i), fn)
        if (err > errmax)
            errmax = err
        end
    end
    return errmax
end

function saveInBucket(err::Float64)
    if err <= 0.5
        return
    end
    ix = 1
    println(err)
    while err > limits[ix]
        println(ix)
        ix += 1
    end
    buckets[ix] += 1
end

function cleanBuckets()
    for i in 1:17
        buckets[i] = 0
    end
end

function saveerrors(fun; path="")
    @showprogress for i in 0:2#(MAX_INT-1)
        x = convert(UInt32, i)
        saveInBucket(check(x, fun))
        #negative values
        x = x + 0x80000000
        saveInBucket(check(x, fun))
    end
    println(fun)
    println(limits)
    println(buckets)
end

for fn in ["sin", "tan", "cos"]
    # for fn in ["csc", "sec", "cot"]
    compile(fn)
    saveerrors(fn; path="/Volumes/TRIG_ERRORS")
    cleanBuckets()
end
# saveerrors(FN)

