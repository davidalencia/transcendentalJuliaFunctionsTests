using Pkg

const zimmermannLib = joinpath(dirname(pathof(ElementalFunctionsTests)), "..", "clibs", "checklib")

#------------------Float64--------------
function ulpdistance(x::Float64, y::BigFloat)::Float64
  return ccall((:ulp_distance, zimmermannLib), Cdouble, (Cdouble, Ref{BigFloat}), x, y)
end

function distance2inf(x::Float64)::Float64
  return ccall((:distance2inf64, zimmermannLib), Cdouble, (Cdouble,), x)
end

function ulperror(foo::Ptr{Cvoid}, mpfr_foo::Ptr{Cvoid}, x::Float64)::Float64
  return ccall((:ulp_error, zimmermannLib), Cdouble,
    (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble), foo, mpfr_foo, x)
end

#------------------Float32--------------
function ulpdistance32(x::Float32, y::BigFloat)::Float64
  return ccall((:ulp_distance32, zimmermannLib), Cdouble, (Cfloat, Ref{BigFloat}), x, y)
end

function distance2inf32(x::Float32)::Float64
  return ccall((:distance2inf32, zimmermannLib), Cdouble, (Cfloat,), x)
end

@inline function ulperror32(foo::Ptr{Cvoid}, mpfr_foo::Ptr{Cvoid}, x::Float32)::Float64
  return ccall((:ulp_error32, zimmermannLib), Cdouble,
    (Ptr{Cvoid}, Ptr{Cvoid}, Cfloat), foo, mpfr_foo, x)
end

#--------------------mpfr----------------
function mpfrfun(strfoo::String)::Ptr{Cvoid}
  @assert strfoo in mpfrfunctions_str "'$strfoo' not in the available mpfr functions"
  return ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), strfoo)
end

@noinline macro ulperror(foo, x)
  @assert string(foo) in mpfrfunctions_str "'$(string(foo))' not in the available mpfr functions"
  quote
    foo_c = @cfunction($foo, Cdouble, (Cdouble,))
    foo_mpfr = ccall((:get_mpfr_fun, zimmermannLib), Ptr{Cvoid}, (Cstring,), string($foo))
    return ccall((:ulp_error, zimmermannLib), Cdouble,
      (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble), foo_c, foo_mpfr, $x)
  end
end