#!/bin/bash
#SBATCH  -J  Heat1D
#SBATCH  -A  eecs
#SBATCH  -p  share
#SBATCH  -N 16      # max number of nodes
#SBATCH  -n 16      # max number of tasks
#SBATCH --constraint=ib
#SBATCH  -o  heat1d.out
#SBATCH  -e  heat1d.err
#SBATCH  --mail-type=BEGIN,END,FAIL
#SBATCH  --mail-user=dengyux@oregonstate.edu
module load openmpi
for s in 1024 8192 65536 524288 4194304 16777216 33554432
do 
    for p in 1 2 4 8 16
    do
        mpic++ -DNUMELEMENTS=$s heat1d_mpi.cpp -o heat1d -lm
	    mpiexec -mca btl self,tcp  -np $p ./heat1d
	    echo " "
    done
done
