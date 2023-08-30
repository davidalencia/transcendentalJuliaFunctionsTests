
export ulpdistance, ulperror, check, get_mpfr_fun

const zimmermannLib = "clibs/checklib.so"
const avaliable_mpfrfunctions = ["cos", "sin", "tan", "cospi", "sinpi", "acos",
  "asin", "atan", "csc", "sec", "cot", "cosh", "sinh", "tanh", "acosh",
  "asinh", "atanh", "exp", "expm1", "exp2", "exp10", "log", "log2",
  "log10", "log1p"]


function ulpdistance(x::Float64, y::BigFloat)::Float64
  return ccall((:ulp_distance, zimmermannLib), Cdouble, (Cdouble, Ref{BigFloat}), x, y)
end

function distance2inf(x::Float64)::Float64
  return ccall((:distance2inf64, zimmermannLib), Cdouble, (Cdouble,), x)
end

function mpfrfun(strfoo::String)::Ptr{Cvoid}
  @assert strfoo in avaliable_mpfrfunctions "'$strfoo' not in the available mpfr functions"
  return ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), strfoo)
end

function ulperror(foo::Ptr{Cvoid}, mpfr_foo::Ptr{Cvoid}, x::Float64)::Float64
  return ccall((:ulp_error, zimmermannLib), Cdouble,
    (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble), foo, mpfr_foo, x)
end

macro ulperror(foo, x)
  @assert string(foo) in avaliable_mpfrfunctions "'$(string(foo))' not in the available mpfr functions"
  @show x
  quote
    foo_c = @cfunction($foo, Cdouble, (Cdouble,))
    foo_mpfr = ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), string($foo))
    return ccall((:ulp_error, zimmermannLib), Cdouble,
      (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble), foo_c, foo_mpfr, $x)
  end
end