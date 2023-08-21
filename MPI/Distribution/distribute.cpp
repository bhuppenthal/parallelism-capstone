#include <mpi.h>

#include "partition.h"

int NumCpus;

// 

void DistributePartitions();
// void DoOneTimeStep();
void GatherResult(int me);

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Status status;

    int me;

    MPI_Comm_Size(MPI_COMM_WORLD, &NumCpus);
    MPI_Comm_Rank(MPI_COMM_WORLD, &me);
}

void DistributePartitions() {
    // Generate the partitions given NUM_PARTITIONS, if using a rectangular partition method
    // Generate_Partitions();

    // Partition the "global array" into local sizes
    // This call will generate a horizontal partition
    Partition_2D_Array(NUM_PARTITIONS, 1);

    // BOSS creates NUM_PARTS 2D arrays


    // other processors create only their own local array to read into

    // send / receive the data, boss can leave in its own

}