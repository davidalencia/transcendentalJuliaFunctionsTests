funs=( "cos" "sin" "tan" )

for i in "${funs[@]}"           
do      
        /opt/homebrew/opt/llvm@16/bin/clang -Wall -shared -O3 -o check.so -lmpfr -lgmp -lm -fopenmp -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -I/opt/homebrew/Cellar/mpfr/4.2.0-p4/include -L/opt/homebrew/Cellar/mpfr/4.2.0-p4/lib check_exhaustive_wrapper_as_argument.c -DSTR=$i
        julia --math-mode=fast test.jl $i
done