#!/bin/bash

SCHEME=rectangular-partition
VERIFY=true
SIMULATE=false

g++ partition.cpp -c -o partition.o
g++ heat.cpp -c -o heat.o
g++ $SCHEME.cpp -DVERIFY_RESULTS=$VERIFY -DPRINT_ALL_TIME_STEPS=$SIMULATE -c -o $SCHEME.o -lm -fopenmp
g++ partition.o heat.o $SCHEME.o -o verify -lm -fopenmp

# Compile times are quick, so we can just remove the object files.
rm -rf *.o

./verify