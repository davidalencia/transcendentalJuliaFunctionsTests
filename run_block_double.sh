funs=( "cos" "sin" "tan" )
ranges=( "-s -pi/4.0 -e pi/4.0" "-s -pi/4.0 -e pi/4.0" "-s -pi/4.0 -e pi/4.0" ) # ⟨ensl-01529804⟩ 9.1.2 Range reduction
# funs=( "cospi" "sinpi")
# ranges=( "-s -1.0 -e 1.0" "-s -1.0 -e 1.0") # ⟨ensl-01529804⟩ 13.1 Range reduction
# funs=( "acos" "asin" "atan" )
# ranges=( "-s -1.0 -e 1.0" "-s -1.0 -e 1.0" "-s 0.0 -e 0x1.fffffffffffffp+1023")
# funs=( "csc" "sec" "cot" )
# ranges=( "-s -pi -e pi --batch 10000000" "-s -pi/2 -e 3*pi/2 --batch 10000000" "-s 0.0 -e pi --batch 10000000")
# funs=( "cosh" "sinh" "tanh" )
# ranges=( "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 1000000" 
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000"
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000")
# funs=( "acosh" "asinh" "atanh" )
# ranges=( "-s 1.0 -e 0x1.fffffffffffffp+1023 --batch 10000000" 
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000"
#  "-s -1 -e +1 --batch 10000000")
# funs=( "exp" "expm1" "exp2" "exp10")
# ranges=( "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000" 
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000"
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000" 
#  "-s -0x1.fffffffffffffp+1023 -e 0x1.fffffffffffffp+1023 --batch 10000000")
# funs=( "log" "log2" "log10" "log1p")
# ranges=( "-s 0.0 -e 0x1.ffffffffffffep+1023 --batch 1000000" 
#  "-s 0.0 -e 0x1.ffffffffffffep+1023 --batch 1000000"
#  "-s 0.0 -e 0x1.ffffffffffffep+1023 --batch 1000000" 
#  "-s -1.0 -e 0x1.ffffffffffffep+1023 --batch 1000000")

for i in "${!funs[@]}"           
do      
        /opt/homebrew/opt/llvm@16/bin/clang -Wall -shared -O3 -o check_double.so -lmpfr -lgmp -lm -I/opt/homebrew/Cellar/mpfr/4.2.0-p12/include -L/opt/homebrew/Cellar/mpfr/4.2.0-p12/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib  check_double.c  -DSTR=${funs[i]}
        julia testdouble.jl ${funs[i]} ${ranges[i]}
        # printf "%s is in %s\n" "${funs[i]}" "${ranges[i]}"
done