#!/bin/bash
for t in 1 2 4 8 12 16 20 24 32
do
  for e in 128 1024 2048 4096 8192 16384 32768 65536 
  do
     g++   openmp-1D.cpp  -DNUMT=$t -DNUME=$e  -o openmp-1D  -lm  -fopenmp
    ./openmp-1D
  done
done