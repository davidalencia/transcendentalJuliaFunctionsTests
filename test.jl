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
    "-I/opt/homebrew/Cellar/mpfr/4.2.0/include",
    "-L/opt/homebrew/Cellar/mpfr/4.2.0/lib",
]
const FILE = "check_exhaustive_wrapper_as_argument.c"
const FN = ARGS[1]


struct UlpError
    ix::Int32
    reinterpret::Float32
    error::Float64
end

function compile(fn)
    cmd = Cmd([COMPILER, FLAGS..., FILE, "-DSTR=$fn"])
    run(cmd)
end

c_f = @cfunction(getfield(Base, Symbol(FN)), Cdouble, (Cdouble,));

function check(x::UInt32)::Float64
    return ccall(
        (:check, "check.so"),
        Cdouble,
        (Cuint, Ptr{Cdouble},),
        x, c_f
    )
end

function getMaxError()
    errmax::Float64 = 0
    err::Float64 = 0
    @showprogress for i in 0:(2_139_095_040-1)
        err = check(convert(UInt32, i))
        if (err > errmax)
            errmax = err
        end
    end
    return errmax
end

function saveerrors(fun; path="")
    open(joinpath(path, "error_$fun.avro"), "w") do f
        @showprogress for i in 0:(MAX_INT-1)
            x = convert(UInt32, i)
            err = UlpError(i, reinterpret(Float32, x), check(x))
            Avro.write(f, err)
            #negative values
            x = x + 0x80000000
            err = err = UlpError(i, reinterpret(Float32, x), check(x))
            Avro.write(f, err)
        end
    end
end

compile(FN)
saveerrors(FN; path="/Volumes/TRIG_ERRORS")
# saveerrors(FN)
