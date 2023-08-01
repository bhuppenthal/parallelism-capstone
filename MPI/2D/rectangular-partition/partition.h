#ifndef PARTITION_H
#define PARTITION_H

// #define NUM_TIME_STEPS  100

#ifndef SIDE
#define SIDE            16
#endif

const int NUME = SIDE*SIDE;

#ifndef NUMCPU
#define NUMCPU            8
#endif

#define NUM_ELEM_PER_THREAD    (NUME/NUMCPU)

#define DEBUG                       true

// Print readable results
#ifndef VERIFY_RESULTS
#define VERIFY_RESULTS              false
#endif

// Print results for the simulation
#ifndef PRINT_ALL_TIME_STEPS
#define PRINT_ALL_TIME_STEPS		false
#endif

// Print CSV formatted performance results
#ifndef CSV
#define CSV                         false
#endif


// Partition information and methods
struct tuple {
    int rows;
    int cols;
};

struct partition {
    int row_start;
    int row_end;
    int col_start;
    int col_end;
};

extern struct partition partitions[NUMCPU];

// Function prototypes.
struct tuple* Generate_Partitions();

void Partition_2D_Array(int Partition_Rows, int Partition_Cols);

void Print_Time_Step(float Temps[2][SIDE][SIDE], int now);

#endif