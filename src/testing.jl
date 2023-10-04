using StatsBase
using Plots
using ProgressMeter
using Random
using Printf
import Base.GC.gc

struct Error
  x::Float64
  err::Float64 
end

struct TestsResults 
  f::Symbol
  bucketsrange::StepRangeLen
  buckets::Vector{Float64}
  maxError::Error
  batchsize::Int
end

const openinfinite = (nextfloat(-Inf), prevfloat(Inf))


buckets(r) = zeros(r.len-1)
hex(x) = @sprintf "%a" x

function add2bucket!(h, r, x)
  if x<r[1] || abs(x) == 0;
    return
  end
  for i in 2:r.len
    if x<r[i] 
      h[i-1]+=1
      return
    end
  end
end

function domain(f::Symbol)::Tuple{Float64,Float64}
  d = Dict(
    :cos => openinfinite,
    :sin => openinfinite,
    :tan => openinfinite,
    :cospi => openinfinite,
    :sinpi => openinfinite,
    :acos => (-1.0, 1.0),
    :asin => (-1.0, 1.0),
    :atan => openinfinite,
    :csc => openinfinite,
    :sec => openinfinite,
    :cot => openinfinite,
    :cosh => openinfinite,
    :sinh => openinfinite,
    :tanh => openinfinite,
    :acosh => (1.0, prevfloat(Inf)),
    :asinh => openinfinite,
    :atanh => (nextfloat(-1.0), prevfloat(1.0)), 
    :exp => (-0x1p20, prevfloat(Inf)), 
    :expm1 => (-0x1p20, prevfloat(Inf)), 
    :exp2 => (-0x1p20, prevfloat(Inf)), 
    :exp10 => (-0x1p20, prevfloat(Inf)), 
    :log => (nextfloat(0.0), prevfloat(Inf)), 
    :log2 => (nextfloat(0.0), prevfloat(Inf)), 
    :log10 => (nextfloat(0.0), prevfloat(Inf)), 
    :log1p => (nextfloat(-1.0), prevfloat(Inf))
  )
  return d[f]
end

function addvalue!(h::Histogram, x::Float64, edges=range(-2, 2, length=81))
  arr = Array{Float64}(undef, 1)
  if abs(x)==0
    return Nothing
  end
  arr[1] = x
  h_tmp = fit(Histogram, arr, edges)
  merge!(h, h_tmp)
  Nothing
end

function init_testRandomBatch(fun::Symbol)
  Random.seed!(20051999)
  fun_c = eval(:@cfunction($fun, Cdouble, (Cdouble,)))
  fun_mpfr = ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), string(fun))
  r = range(-2, 2, length=41)
  h = buckets(r)
  maxerror = Error(0.0, 0.0)
  return (fun_c, fun_mpfr, h, r, maxerror)
end

@noinline function testRandomBatch(fun::Symbol, rng::AbstractRNG; batchSize=100_000, io = IOBuffer())
  println(io, fun)
  fun_c, fun_mpfr, h, r, maxerror = init_testRandomBatch(fun)
  @showprogress 0.5 string(fun) for i in 1:batchSize
    rn = rand(rng)
    err = ulperror(fun_c, fun_mpfr, rn)
    add2bucket!(h, r, err)
    if err>maxerror.err
      maxerror = Error(rn, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  return TestsResults(fun,r,h, maxerror, batchSize)
end

function save_testRandomBatch(fun::Symbol, rng::AbstractRNG; file ="results/results.txt", imagesDir="./", kwargs...)
  open(file, "a+") do file
    res = testRandomBatch(fun, rng; io=file, kwargs...)
    savefig(bar(res.bucketsrange, res.buckets, label="$fun-$(res.batchsize)"), imagesDir*"$fun-$batchSize.png")
  end
  gc()
end