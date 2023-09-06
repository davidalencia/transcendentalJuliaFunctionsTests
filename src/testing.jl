using StatsBase
using Plots
using ProgressMeter
using Random
import Base.GC.gc

struct Error
  x::Float64
  err::Float64 
end

const openinfinite = (nextfloat(-Inf), prevfloat(Inf))

function add2bucket!(h, r, x)
  if x<r[1]
    return
  end
  for i in 2:r.len
    if x<r[i] 
      h[i-1]+=1
      return
    end
  end
end


buckets(r) = zeros(r.len-1)
hex(x) = @sprintf "%a" x

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
  h = Histogram(range(-2, 2, length=81))
  maxerror = Error(0.0, 0.0)
  return (fun_c, fun_mpfr, h, maxerror)
end

@noinline function testRandomBatch(fun; batchSize=100_000, rand=() -> γsectionCC(domain(Symbol(fun))...), io = stdout::IO, imagesDir="results/")
  println(io, fun)
  fun_c, fun_mpfr, h, maxerror = init_testRandomBatch(fun)
  @showprogress 0.5 string(fun) for i in 1:batchSize
    r = rand()
    err = ulperror(fun_c, fun_mpfr, r)
    addvalue!(h, err)
    if err>maxerror.err
      maxerror = Error(r, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  savefig(plot(h, label="$fun-$batchSize-$(hex.(domain(Symbol(fun))))"), imagesDir*"$fun-$batchSize.png")
end

function save_testRandomBatch(fun; file ="results/results.txt", kwargs...)
  open(file, "a+") do file
    testRandomBatch(fun; io=file, kwargs...)
  end
  gc()
end


function testRandomBatch_γ(fun::Function, domain::Tuple{Float64, Float64}; batchSize::Int=100_000, io::IO = stdout::IO, imagesDir::String="results/")
  println(io, fun)
  fun_c, fun_mpfr, h, maxerror = init_testRandomBatch(fun)
  @showprogress 0.5 string(fun) for i in 1:batchSize
    r = γsectionCC(domain...)
    err = ulperror(fun_c, fun_mpfr, r)
    addvalue!(h, err)
    if err>maxerror.err
      maxerror = Error(r, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  savefig(plot(h, label="$fun-$batchSize-$(hex.(domain(Symbol(fun))))"), imagesDir*"$fun-$batchSize.png")
end

function testRandomBatch_DSR(fun::Symbol, domain::Tuple{Float64, Float64}; batchSize::Int=100_000, io::IO = stdout::IO, imagesDir::String="results/")
  println(io, fun)
  fun_c, fun_mpfr, h, maxerror = init_testRandomBatch(fun)
  @showprogress for i in 1:batchSize
    r = DSR(domain...)
    err = ulperror(fun_c, fun_mpfr, r)
    addvalue!(h, err)
    if err>maxerror.err
     maxerror = Error(r, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  savefig(plot(h, label="$fun-$batchSize-$(hex.(domain(Symbol(fun))))"), imagesDir*"$fun-$batchSize.png")
end

@noinline function testRandomBatch_(fun; batchSize=100_000, rand=() -> γsectionCC(domain(Symbol(fun))...), io = stdout::IO, imagesDir="results/")
  println(io, fun)
  fun_c, fun_mpfr, _, maxerror = init_testRandomBatch(fun)
  r = range(-2, 2, length=81)
  h = buckets(r)
  @showprogress 0.5 string(fun) for i in 1:batchSize
    rn = rand()
    err = ulperror(fun_c, fun_mpfr, rn)
    add2bucket!(h, r, err)
    if err>maxerror.err
      maxerror = Error(rn, err)
    end
  end
  println(io, "maxerror = $(maxerror.err) at $(hex(maxerror.x))")
  println(io, h)
  # savefig(bar(r, h, label="$fun-$batchSize"), imagesDir*"$fun-$batchSize.png")
end
