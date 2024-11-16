#!/bin/bash
export OMP_NUM_THREADS=16

echo "Setting up main..."
GSL=$(brew --prefix gsl)
IFLAGS="-I$GSL/include"
LFLAGS="-Wl,-rpath,$GSL/lib -L$GSL/lib -lm -lgsl"

# LFLAGS="-Wl,-rpath,$GSL/lib -L$GSL/lib -lm -lomp -lgsl"
# CC="mpicxx -Xpreprocessor -fopenmp -O3 -std=c++17"
CC="clang++ -O3 -std=c++17"
CC="clang++ -O3 -std=c++17 -I/opt/homebrew/Cellar/gsl/2.7.1/include"
$CC $IFLAGS -c main.cc

#$CC $IFLAGS -c module1.cc
# $CC $IFLAGS -c module2.cc
$CC $LFLAGS -o main main.cc #Other scripts
echo   # Blank line
echo "Time for execution.."
time ./main # | tee main.out

echo "Starting process..."

rm main.o main
