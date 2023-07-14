#ifndef PARTITION_H
#define PARTITION_H

#define NUM_TIME_STEPS  100                      // number of time steps the simultation runs through

#ifndef SIDE
#define SIDE            8
#endif

const int NUME = SIDE*SIDE;

#ifndef NUMT
#define NUMT            4                       // number of threads to use
#endif

#define NUM_ELEM_PER_THREAD    (NUME/NUMT)      // number of elements in each thread

#ifndef PRINT_ALL_TIME_STEPS
#define PRINT_ALL_TIME_STEPS		true         // set to true to allow all time steps to print
#endif

// Partition information and methods
struct partition {
    int row_start;
    int row_end;
    int col_start;
    int col_end;
};

extern struct partition partitions[NUMT];

// Function prototypes.
void Partition_2D_Array(int Partition_Rows, int Partition_Cols);

#endif