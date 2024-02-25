using StatsBase
using Plots
using ProgressMeter
using Random
using Printf
import Base.GC.gc
import Base.@kwdef


abstract type AbstractError end

@kwdef struct Error <: AbstractError
  x::Float64
  err::Float64 
end

@kwdef struct Error32 <: AbstractError
  x::Float32
  err::Float64 
end

@kwdef struct TestsResults
  f::Symbol
  bucketsrange::StepRangeLen
  buckets::Vector{Float64}
  maxError::AbstractError
  batchsize::Int
  mean::Float64
  variance::Float64
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
  fun_c = eval(:@cfunction($fun, Cdouble, (Cdouble,)))
  fun_mpfr = ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), string(fun))
  r = range(-2, 2, length=41)
  h = buckets(r)
  maxerror = Error(0.0, 0.0)
  return (fun_c, fun_mpfr, h, r, maxerror)
end

function testRandomBatch(fun::Symbol, rng::AbstractRNG; batchSize=100_000, io = IOBuffer(),  seed=20051999)
  println(io, fun)
  Random.seed!(seed)
  sum=0    # To store sum of stream 
  sumsq=0  # To store sum of square of stream 
  n=0   
  fun_c, fun_mpfr, h, r, maxerror = init_testRandomBatch(fun)
  @showprogress 0.5 string(fun) for _ in 1:batchSize
    rn = rand(rng)
    err = ulperror(fun_c, fun_mpfr, rn)
    n+=1
    abserr = abs(err)
    sum+=abserr
    sumsq+=(abserr*abserr) 
    add2bucket!(h, r, err)
    if abserr>abs(maxerror.err)
      maxerror = Error(rn, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  mean = sum/n 
  var = (sumsq/n) - (mean*mean) 
  return TestsResults(fun,r,h, maxerror, batchSize, mean, var)
end

function save_testRandomBatch(fun::Symbol, rng::AbstractRNG; file ="results/results.txt", imagesDir="./", kwargs...)
  open(file, "a+") do file
    res = testRandomBatch(fun, rng; io=file, kwargs...)
    savefig(bar(res.bucketsrange, res.buckets, label="$fun-$(res.batchsize)"), imagesDir*"$fun-$batchSize.png")
  end
  gc()
end