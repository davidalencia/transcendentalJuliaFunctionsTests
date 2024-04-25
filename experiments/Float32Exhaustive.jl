using Pkg
Pkg.activate("..")
using ElementalFunctionsTests

using ProgressMeter
using JSON
using Printf

JSON.lower(e::Error32) =Dict(:x=>@sprintf("%a", e.x), :err=>e.err)

const r = range(-2, 2, length=41)
errors = []

for f in mpfrfunctions
    x = nextfloat(-Inf32)
    maxerror = Error32(0.0, 0.0)
    h = buckets(r)
    
    f_c = eval(:(@cfunction($f, Cfloat, (Cfloat,))))
    f_mpfr = mpfrfun(string(f))
    dom = domain(Symbol(f))
    prog = Progress(4278190078, showspeed=true, dt=2.0, desc=string(f))
    sum=0    # To store sum of stream 
    sumsq=0  # To store sum of square of stream 
    n=0   

    
    while(!isinf(x))
        next!(prog)
        if x < dom[1]
            x = nextfloat(x)
            continue
        elseif x > dom[2]
            break;
        end
        
        err = ulperror32(f_c, f_mpfr, x)
        add2bucket!(h, r, err)
        abserr = abs(err)
        n+=1
        sum+=abserr
        sumsq+=(abserr*abserr) 
        
        if abs(err)>abs(maxerror.err)
          maxerror = Error32(x, err)
        end
        x = nextfloat(x)
    end

    mean = sum/n 
    var = (sumsq/n) - (mean*mean) 
    
    finish!(prog)
      
    push!(errors, TestsResults(Symbol(f),r,h, maxerror, 0, mean, var))
    open("Float32Errors.json","w") do f
      write(f, JSON.json(errors))
    end
end
