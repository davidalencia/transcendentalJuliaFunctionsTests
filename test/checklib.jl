using Test
using PyCall

if Sys.iswindows()
    const libmpfr = "libmpfr-6.dll"
elseif Sys.isapple()
    const libmpfr = "@rpath/libmpfr.6.dylib"
else
    const libmpfr = "libmpfr.so.6"
end

hex(x) = py"float.hex"(x)
nextinbig(x, step=1) = BigFloat(hex(nextfloat(x, step)))
previnbig(x, step=1) = BigFloat(hex(prevfloat(x, step)))

@testset verbose=true "ulp_distance" begin
    ulpdistance(x, y) = ccall((:ulp_distance, "checklib.so"), Cdouble, (Cdouble, Ref{BigFloat}), x, y)

    @testset "equality" begin 
        @test ulpdistance(0.0, BigFloat("0.0")) == 0
        @test ulpdistance(0x2p-1023, BigFloat("0x2p-1023")) ≈ 0.0
        @test ulpdistance(0x2p-1025, BigFloat("0x2p-1025")) ≈ 0.0
    end

    @testset "directions" begin 
        @test sign(ulpdistance(1.0, nextinbig(1.0))) == 1.0
        @test sign(ulpdistance(1.0, previnbig(1.0))) == -1.0
    end

    @testset "normal power of 2" begin 
        @test ulpdistance(1.0, nextinbig(1.0)) ≈ 1.0
        @test ulpdistance(1.0, previnbig(1.0)) ≈ -0.5
        @test ulpdistance(2.0, nextinbig(2.0)) ≈ 1.0
        @test ulpdistance(prevfloat(2.0), nextinbig(2.0)) ≈ 3.0
        @test ulpdistance(0x2p-1023, nextinbig(0x2p-1023)) ≈ 1.0
        @test ulpdistance(0x2p-1023, previnbig(0x2p-1023)) ≈ -1.0

    end

    @testset "normal not power of 2" begin 
        @test ulpdistance(1.5, nextinbig(1.5)) ≈ 1.0
        @test ulpdistance(1.5, previnbig(1.5)) ≈ -1.0
        @test ulpdistance(5.0, nextinbig(5.0)) ≈ 1.0
        @test ulpdistance(5.0, nextinbig(5.0, 100)) ≈ 100.0
    end

    @testset "subnormals" begin 
        @test ulpdistance(0x2p-1024, nextinbig(0x2p-1024)) ≈ 1.0
        @test ulpdistance(0x2p-1028, nextinbig(0x2p-1028)) ≈ 1.0
        @test ulpdistance(0x2p-1075, nextinbig(0x2p-1075)) ≈ 1.0
        @test ulpdistance(0x2p-1028, nextinbig(0x2p-1028, 100)) ≈ 100.0
    end
    
    @testset "not full ulps" begin 
        @test ulpdistance(1.0, BigFloat(1.0)+BigFloat(hex(eps(1.0)/2))) ≈ 0.5
        @test ulpdistance(1.0, BigFloat(1.0)+BigFloat(hex(eps(1.0)/4))) ≈ 0.25
        @test ulpdistance(1.0, BigFloat(1.0)+BigFloat(hex(eps(1.0)/8))) ≈ 0.125
    end
end

@testset "distance2inf" begin
    distance2inf(x::Float64) = ccall((:distance2inf64, "checklib.so"), Cdouble, (Cdouble,), x)
    
    @test distance2inf(Inf) == 0.0
    @test distance2inf(-Inf) == 0.0
    @test distance2inf(prevfloat(Inf)) == 1.0
    @test distance2inf(nextfloat(-Inf)) == 1.0
    @test distance2inf(prevfloat(Inf, 100)) == 100.0
    @test distance2inf(nextfloat(-Inf, 100)) == 100.0
end


function ulp_error_call(foo::Ptr{Cvoid}, mpfr_foo::Ptr{Cvoid}, x::Float64=100.0)
    return ccall((:ulp_error, "checklib.so"), Cdouble, (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble), foo, mpfr_foo, x)
end

@testset "ulp_error" begin 

    I_c = @cfunction(Cdouble, (Cdouble,)) do x
        return x
    end
    oneulp_c = @cfunction(Cint, (Ref{BigFloat}, Ref{BigFloat}, Cint)) do rop, op, rnd
        y = nextfloat(Float64(op))
        ret = ccall((:mpfr_set_d, libmpfr), Cint, (Ref{BigFloat}, Cdouble, Cint), rop, y, 0)
        if ret!=0 
            throw(ErrorException("something went wrong with mpfr_set_d"))
        end
        return Int32(0) 
    end;
    tenulp_c = @cfunction(Cint, (Ref{BigFloat}, Ref{BigFloat}, Cint)) do rop, op, rnd
        y = nextfloat(Float64(op), 10)
        ret = ccall((:mpfr_set_d, libmpfr), Cint, (Ref{BigFloat}, Cdouble, Cint), rop, y, 0)
        if ret!=0 
            throw(ErrorException("something went wrong with mpfr_set_d"))
        end
        return Int32(0) 
    end;
    inf_onestep_c = @cfunction(Cint, (Ref{BigFloat}, Ref{BigFloat}, Cint)) do rop, op, rnd
        y = prevfloat(Inf)
        ret = ccall((:mpfr_set_d, libmpfr), Cint, (Ref{BigFloat}, Cdouble, Cint), rop, y, 1) ## 1  = RNDZ
        println("Julia] ret=",ret)
        # if ret>0
        #     throw(ErrorException("something went wrong with mpfr_set_d"))
        # end
        return Int32(0) 
    end;
    inf_tenstep_c = @cfunction(Cint, (Ref{BigFloat}, Ref{BigFloat}, Cint)) do rop, op, rnd
        y = prevfloat(Inf, 10)
        ret = ccall((:mpfr_set_d, libmpfr), Cint, (Ref{BigFloat}, Cdouble, Cint), rop, y, 0)
        println("Julia] ret=",ret)
        if ret!=0 
            throw(ErrorException("something went wrong with mpfr_set_d"))
        end
        return Int32(0) 
    end;


    @test ulp_error_call(I_c, oneulp_c) ≈ 1.0
    @test ulp_error_call(I_c, tenulp_c) ≈ 10.0
    @test ulp_error_call(I_c, inf_tenstep_c, Inf) ≈ 10.0
    @test ulp_error_call(I_c, inf_onestep_c, Inf) ≈ 1.0

end
