using ProgressMeter
using Plots, StatsBase
using Distributions
using PyCall
import Random
using ArgParse

abstract type ParsedFloat <: AbstractFloat end

function ArgParse.parse_item(::Type{Float64}, x::AbstractString)
    return Float64(eval(Meta.parse(x)))
end

s = ArgParseSettings()
@add_arg_table s begin
  "fun"
      required = true
  "--rng_start", "-s"
      required = true
      arg_type = Float64
  "--rng_end", "-e"
      required = true
      arg_type = Float64
  "--batch"
      arg_type = Int
      default = 1000000
end

args = parse_args(ARGS, s)

const FN = get(args, "fun", "exp")
const RNG_START = get(args, "rng_start", -1.0)
const RNG_END = get(args, "rng_end", 1.0)
const BATCH = get(args, "batch", 1000)

const positiverange = 0x0:0x7ff0000000000000
const negativerange = 0x8000000000000000:0xfff0000000000000
const edges = -1.5:0.1:1.5

const c_fun = @cfunction(getfield(Main, Symbol(FN)), Cdouble, (Cdouble,))
const rng = (RNG_START, RNG_END)

function hex(x)
    return py"float.hex"(x)
end

function check(x::UInt64)::Float64
    return ccall(
        (:check, "check_double.so"),
        Cdouble,
        (Culong, Ptr{Cdouble},),
        x, c_fun
    )
end

function _getRandomFloat(a,b)
    a = reinterpret(UInt64, a)
    b = reinterpret(UInt64, b)
    @assert a<b "a must be less than b in UInt form"
    c = floor(rand()*(b-a)+a)
    c = UInt64(c)
    return reinterpret(Float64, c)
end

function getRandomFloat(a, b) #not a random number a random float
    @assert a<b "a must be less than b"
    if a<0 && b<=0
        return _getRandomFloat(b,a) # when using UInt -0.0 is less than -Inf
    elseif a>=0 && b>0
        return _getRandomFloat(a,b)
    else
       a = _getRandomFloat(-0.0, a) # when using UInt -0.0 is less than -Inf
       b = _getRandomFloat(+0.0 ,b)
       return rand((a,b))
    end 

end

function montecarlo(r, poolsize::Int)
    max = (-1.0, -1, 0.0)
    Random.seed!(187372311)
    i = 0
    h = fit(Histogram, [], edges)
    results = []

    @showprogress for i in 1:poolsize
        x = getRandomFloat(r[1], r[2])
        num = reinterpret(UInt64, x)
        err = check(num)
        if(abs(err)>0)
            push!(results, err)
        end
        i += 1
        if (i>100)
            h_tmp = fit(Histogram, results , edges)
            results = []
            merge!(h, h_tmp)
        end
        if(abs(err)>abs(max[3]))
            max=(x, num, err)
        end
    end
    # println(FN)
    # println(max)
    savefig(plot(h, label="$FN-$poolsize"), "results/$FN-$poolsize.png")
    file = open("results/results.txt", "a+")

    write(file, "$FN-$poolsize\n")
    write(file, "max error $(max[3]) on $(hex(max[1]))\n")
    write(file, "$(h.weights) \n")
    write(file, "-------------------------------------------------------\n")

    close(file)

end

montecarlo(rng , BATCH) 

