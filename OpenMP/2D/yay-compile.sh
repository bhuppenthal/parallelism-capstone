#!/bin/bash

g++ partition.cpp -c -o partition.o
g++ horizontal-partition.cpp -c -o horizontal-partition.o -lm -fopenmp
g++ partition.o horizontal-partition.o -o horizontal-verify -lm -fopenmp