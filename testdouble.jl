using ProgressMeter
import Random

# +0_64 = 0
# +Inf_64 = 9218868437227405312

const FN = ARGS[1]
const positiverange = 0x0:0x7ff0000000000000
const negativerange = 0x8000000000000000:0xfff0000000000000

const c_fun = @cfunction(getfield(Main, Symbol(FN)), Cdouble, (Cdouble,))

function check(x::UInt64)::Float64
    return ccall(
        (:check, "check_double.so"),
        Cdouble,
        (Culong, Ptr{Cdouble},),
        x, c_fun
    )
end

function montecarlo(r::UnitRange{UInt}, poolsize::Int)
    # println(r.start)
    # println(r.stop)
    # for i in 
    max = (-1, -1, 0.0)
    Random.seed!(187372311)

    @showprogress for i in 1:poolsize
        x = rand(r)
        err = check(x)
        num = reinterpret(Float64, x)
        # println("$(num) : $(err)")
        if(err>max[3])
            max=(x, num, err)
            println(max)
        end
    end
    println(max)
    # println(reinterpret(Float64, rand(r, 100)))
end

montecarlo(positiverange, 100_000_000_000) 
# 0x0-0x7ff0000000000000, 0x7ff0000000000001-0x7fffffffffffffff, 0x8000000000000000-0xfff0000000000000, 0xfff0000000000001-0xffffffffffffffff
# +0     +Inf                   NaN            NaN                     -0                   -Inf            NaN                     NaN