#!/bin/bash

# Uncomment lines 4-13 to check performance values of OpenMP 2D horizontal parition
# for t in 1 2 4 8 12 16 20 24 32
# do
#   for s in 32 96 256 512 1024 2048 4096
#   do
#     g++ partition.cpp -DNUMT=$t -DSIDE=$s -c -o partition.o
#     g++ horizontal-partition.cpp -DNUMT=$t -DSIDE=$s -c -o horizontal-partition.o -lm -fopenmp
#     g++ partition.o horizontal-partition.o -DNUMT=$t -DSIDE=$s -o horizontal-verify -lm -fopenmp
#     ./horizontal-verify
#   done
# done

#Uncomment lines 16-25 to check performance values of OpenMP 2D vertical paritions - column major order
# for t in 1 2 4 8 12 16 20 24 32
# do
#   for s in 32 96 256 512 1024 2048
#   do
#     g++ partition.cpp -DNUMT=$t -DSIDE=$s -c -o partition.o
#     g++ vertical-partition-col-traversal.cpp -DNUMT=$t -DSIDE=$s -c -o vertical-partition-col-traversal.o -lm -fopenmp
#     g++ partition.o vertical-partition-col-traversal.o -DNUMT=$t -DSIDE=$s -o vertical-col-verify -lm -fopenmp
#     ./vertical-col-verify
#   done
# done

#Uncomment lines 28-37 to check performance values of OpenMP 2D vertical paritions - row major order
# for t in 1 2 4 8 12 16 20 24 32
# do
#   for s in 32 96 256 512 1024 2048
#   do
#     g++ partition.cpp -DNUMT=$t -DSIDE=$s -c -o partition.o
#     g++ vertical-partition-row-traversal.cpp -DNUMT=$t -DSIDE=$s -c -o vertical-partition-row-traversal.o -lm -fopenmp
#     g++ partition.o vertical-partition-row-traversal.o -DNUMT=$t -DSIDE=$s -o vertical-row-verify -lm -fopenmp
#     ./vertical-row-verify
#   done
# done