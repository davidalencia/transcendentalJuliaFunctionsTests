using Distributions
using Printf

export randomfloat, γsectionCC


"""
TODO: modificar γ para que si a<0 & b>0 & abs(a)>abs(b) entonces has lo que hice
sino toma ulp(b)
tldr: toma el ulp mas grande dentro del intervalo.

"""
@inline function γ(a, b)
  a_abs = abs(a)
  b_abs = abs(b)
  a_ulp = nextfloat(a_abs) - a_abs
  b_ulp = b_abs - prevfloat(b_abs)
  return max(a_ulp, b_ulp)
end

"""
ceilint(a,b,g)
Compute db/g-a/ge correctly using Dekker’s algorithm for an exact summation. 
"""
@inline function ceilint(a, b, g)::Int64
  s = b / g - a / g
  if abs(a) <= abs(b)
    ε = -a / g - (s - b / g)
  else
    ε = b / g - (s + a / g)
  end
  si = ceil(Int, s)
  return (s != si) ? si : si + Int(ε > 0)
end

"""
F. Goualard γsectionCC(a,b)
γ section algorithm for close close intervals
Draw a float from an interval [a,b] uniformly at random.
"""
function γsectionCC(a, b)
  g = γ(a, b)
  hi = ceilint(a, b, g)
  k = rand(DiscreteUniform(0, hi))
  if abs(a) <= abs(b)
    return (k == hi) ? a : b - k * g
  else
    return (k == hi) ? b : a + k * g
  end
end

"""
Same Sign Random
"""
function SSR(a::Float64, b::Float64)::Float64
  a = reinterpret(UInt64, a)
  b = reinterpret(UInt64, b)
  a, b = a < b ? (a, b) : (b, a)
  c = rand(DiscreteUniform(a, b))
  c = UInt64(c)
  return reinterpret(Float64, c)
end

function dis(a::Float64, b::Float64)::UInt64
  return reinterpret(UInt64, b) - reinterpret(UInt64, a)
end

function DSR(a::Float64, b::Float64)::Float64
  @assert a < b "a must be less than b"
  if a == 0.0
    return SSR(+0.0, b)
  elseif b == 0.0
    return SSR(-0.0, a)
  elseif sign(a) == sign(b)
    return SSR(a, b)
  else
    l = dis(-0.0, a)
    r = dis(+0.0, b)
    k = rand(DiscreteUniform(-Int64(l), Int64(r)))
    if k < 0
      return reinterpret(Float64, reinterpret(UInt64, -0.0) + UInt64(-k))
    else
      return reinterpret(Float64, UInt64(k))
    end
  end
end


