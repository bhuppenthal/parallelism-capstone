#!/bin/bash

#define which file to test and output type
SCHEME=rectangular-partition
VERIFY=true

NUMT=32
SIDE=8

uptime

g++ heat.cpp -c -o heat.o
g++ partition.cpp -DNUMT=$NUMT -DSIDE=$SIDE -c -o partition.o
g++ $SCHEME.cpp -DVERIFY_RESULTS=$VERIFY -DNUMT=$NUMT -DSIDE=$SIDE -c -o $SCHEME.o -lm -fopenmp
g++ partition.o heat.o $SCHEME.o -o verify -lm -fopenmp

./verify

rm -rf *.o
rm verify