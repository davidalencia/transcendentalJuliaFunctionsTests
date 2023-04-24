funs=( "cos" "sin" "tan" )

for i in $funs           
do
        julia test.jl $i
done