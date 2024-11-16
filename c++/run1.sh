# #!/bin/bash

# function compile() {
#     compiler="/opt/homebrew/bin/mpicxx -O3 -std=c++17"
#     $compiler -c $1.cc
#     $compiler -o $1 $1.o -lm -lgsl
#     rm $1.o
# }

# function run() {
#     compile $1
#     time ./$1
#     rm $1
# }

# run solve_levels

#!/bin/bash
export OMP_NUM_THREADS=1

echo "Setting up main..."
GSL=$(brew --prefix gsl)
IFLAGS="-I$GSL/include"
LFLAGS="-Wl,-rpath,$GSL/lib -L$GSL/lib -lm -lgsl"

# LFLAGS="-Wl,-rpath,$GSL/lib -L$GSL/lib -lm -lomp -lgsl"
# CC="mpicxx -Xpreprocessor -fopenmp -O3 -std=c++17"
CC="clang++ -O3 -std=c++17"
CC="clang++ -O3 -std=c++17 -I/opt/homebrew/Cellar/gsl/2.7.1/include"
$CC $IFLAGS -c solve_levels.cc

#$CC $IFLAGS -c module1.cc
# $CC $IFLAGS -c module2.cc
$CC $LFLAGS -o solve_levels solve_levels.cc #main main.cc #Other scripts
echo   # Blank line
echo "Time for execution.."
time ./main # | tee main.out

echo "Starting process..."

rm solve_levels
