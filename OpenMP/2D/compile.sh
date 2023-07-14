#!/bin/bash

g++ partition.cpp -c -o partition.o
g++ heat.cpp -c -o heat.o
g++ horizontal-partition.cpp -c -o horizontal-partition.o -lm -fopenmp
g++ partition.o heat.o horizontal-partition.o -o verify -lm -fopenmp