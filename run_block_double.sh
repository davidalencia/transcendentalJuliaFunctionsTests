# funs=( "cos" "sin" "tan" )
# ranges=( "-pi+0.0 pi+0.0" "0.0 2pi" "-pi/2 pi/2")
# funs=( "csc" "sec" "cot" )
# ranges=( "-pi+0.0 pi+0.0" "-pi/2 3*pi/2" "0.0 pi+0.0")
# funs=( "acos" "asin" "atan" )
# ranges=( "-1.0 1.0" "-1.0 1.0" "0.0 0x1.ffffffffffffep+1023")
funs=( "exp" )
ranges=( "-0x1.ffffffffffffep+1023 0x1.ffffffffffffep+1023")

for i in "${!funs[@]}"           
do      
        /opt/homebrew/opt/llvm@16/bin/clang -Wall -shared -O3 -o check_double.so -lmpfr -lgmp -lm -I/opt/homebrew/Cellar/mpfr/4.2.0-p9/include -L/opt/homebrew/Cellar/mpfr/4.2.0-p9/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib  check_double.c  -DSTR=${funs[i]}
        julia testdouble.jl ${funs[i]} ${ranges[i]}
        # printf "%s is in %s\n" "${funs[i]}" "${ranges[i]}"
done