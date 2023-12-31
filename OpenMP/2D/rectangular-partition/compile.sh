#!/bin/bash

#define which file to test and output type
SCHEME=rectangular-partition
DEBUG=false
VERIFY=false
SIMULATE=false
CSV=true
OUTPUT=rectangular-test.csv

uptime

g++ heat.cpp -c -o heat.o

for NUMT in 1 2 4 8 12 16 20 24 32
do
    for SIDE in 32 96 256 512 1024 2048 4096
    do
        g++ partition.cpp -DDEBUG=$DEBUG -DNUMT=$NUMT -DSIDE=$SIDE -c -o partition.o
        g++ $SCHEME.cpp -DDEBUG=$DEBUG -DVERIFY_RESULTS=$VERIFY -DPRINT_ALL_TIME_STEPS=$SIMULATE -DCSV=$CSV -DNUMT=$NUMT -DSIDE=$SIDE -c -o $SCHEME.o -lm -fopenmp
        g++ partition.o heat.o $SCHEME.o -o verify -lm -fopenmp

        rm -rf partition.o $SCHEME.o

        ./verify
    done
done

rm -rf *.o
rm verify