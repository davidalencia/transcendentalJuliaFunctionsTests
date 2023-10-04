import Base.rand
import Random.AbstractRNG

export Floating64Distribution, FloatingDistribution
export γSection64, γSection

#------------------------------------------------------------------------------
#----------------------------------structs-------------------------------------
#------------------------------------------------------------------------------

struct Floating64Distribution <: AbstractRNG
  a::Float64
  b::Float64
  Floating64Distribution(a::Float64, b::Float64) =
    a < b ? new(a, b) : error("a must be less than b")
end

const FloatingDistribution = Floating64Distribution

struct γSection64 <: AbstractRNG
  a::Float64
  b::Float64
  g::Float64
  hi::Int64
  function γSection64(a::Float64, b::Float64)
    @assert a < b "a must be less than b"
    g = γ(a, b)
    return new(a, b, g, ceilint(a, b, g))
  end
end

const γSection = γSection64


#------------------------------------------------------------------------------
#-------------------------------Util functions---------------------------------
#------------------------------------------------------------------------------


"""
tldr: toma el ulp mas grande dentro del intervalo.
"""
function γ(a, b)
  a_ulp = nextfloat(a) - a
  b_ulp = b - prevfloat(b)
  return max(a_ulp, b_ulp)
end

"""
ceilint(a,b,g)
Compute b/g-a/g correctly using Dekker’s algorithm for an exact summation. 
"""
function ceilint(a, b, g)::Int64
  s = b / g - a / g
  if abs(a) <= abs(b)
    ε = -a / g - (s - b / g)
  else
    ε = b / g - (s + a / g)
  end
  si = ceil(Int, s)
  return (s != si) ? si : si + Int(ε > 0)
end

function samesignRandom(a::UInt64, b::UInt64)::Float64
  return reinterpret(Float64, rand(a:b))
end


#------------------------------------------------------------------------------
#------------------------------------rand--------------------------------------
#------------------------------------------------------------------------------

@inline function rand(rng::Floating64Distribution)::Float64
  a_uint = reinterpret(UInt64, rng.a)
  b_uint = reinterpret(UInt64, rng.b)
  zp_uint = reinterpret(UInt64, +0.0)
  zn_uint = reinterpret(UInt64, -0.0)
  if rng.a == 0.0
    return samesignRandom(zp_uint, b_uint)
  elseif rng.b == 0.0
    return samesignRandom(zn_uint, a_uint)
  elseif sign(rng.a) == sign(rng.b)
    a_uint, b_uint = a_uint < b_uint ? (a_uint, b_uint) : (b_uint, a_uint)
    return samesignRandom(a_uint, b_uint)
  end
  l = a_uint - zn_uint
  k = rand(-Int64(l):Int64(b_uint))
  if k < 0
    return reinterpret(Float64, zn_uint + UInt64(-k))
  else
    return reinterpret(Float64, UInt64(k))
  end
end


@inline function rand(rng::γSection64)::Float64
  k = rand(0:rng.hi)
  kg = Float64(k) * rng.g
  top = k == rng.hi
  cmp_abs = abs(rng.a) <= abs(rng.b)
  if isinf(kg)
    kg = ((k / 2) * (rng.g / 2)) * 2
  end
  if top && cmp_abs
    return rng.a
  elseif top
    return rng.b
  elseif cmp_abs
    return rng.b - kg
  else
    return rng.a + kg
  end
end
