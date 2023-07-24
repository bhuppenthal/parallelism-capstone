#!/bin/bash
#SBATCH  -J  Heat2D
#SBATCH  -A  eecs
#SBATCH  -p  share
#SBATCH  -N 16      # max number of nodes
#SBATCH  -n 16      # max number of tasks
#SBATCH --constraint=ib
#SBATCH  -o  heat2d.out
#SBATCH  -e  heat2d.err
#SBATCH  --mail-type=END,FAIL
#SBATCH  --mail-user=immermam@oregonstate.edu
module load openmpi
mpic++ heat2d.cpp -o heat2d -lm
for p in 1 4 9 16
do
	mpiexec -mca btl self,tcp -np $p ./heat2d
	echo " "
done
