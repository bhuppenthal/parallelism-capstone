#!/bin/bash
#SBATCH  -J  Heat2D
#SBATCH  -A  eecs
#SBATCH  -p  share
#SBATCH  -N 16      # max number of nodes
#SBATCH  -n 16      # max number of tasks
#SBATCH --constraint=ib
#SBATCH  -o  heat2d.out
#SBATCH  -e  heat2d.err
#SBATCH  --mail-type=BEGIN,END,FAIL
#SBATCH  --mail-user=getekj@oregonstate.edu
module load openmpi
for s in 32 96 256 512 1024 2048 4096 6144
do 
    for p in 1 2 4 8 16
    do
        mpic++ -DGRID_SIZE=$s -DNUMCPU=$p -DSIDE=$s heat.cpp partition.cpp rectangular-partition.cpp -o rectangular-heat2d -lm
	    mpiexec -mca btl self,tcp  -np $p ./rectangular-heat2d
	    echo " "
    done
done