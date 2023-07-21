/*
    Heat Diffusion Equation Simulation - 2D Using OpenMP

    Naive implementation using horizontal partitions.
*/

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

#include "partition.h"
#include "heat.h"


float   Temps[2][SIDE][SIDE];

int     Now;                                    // which array is the "current values" = 0 or 1
int     Next;                                   // which array is being filled = 1 or 0

void    DoAllWork(int);

int main(void) {

    Now = 0;
    Next = 1;

    // Setting number of threads that will be used
    omp_set_num_threads(NUMT);

    // Setting all initial temperatures to 0, except the middle value which is set to 100
    for (int i = 0; i < SIDE; i++)
        for (int j = 0; j < SIDE; j++)
            Temps[Now][i][j] = (float) 0;
    
    Temps[Now][SIDE/2][SIDE/2] = (float) 100;

    if (VERIFY_RESULTS) {
        for (int i = 0; i < SIDE; i++) {
            for (int j = 0; j < SIDE; j++)
                printf(" %.2f ", Temps[Now][i][j]);
            printf("\n");
        }
    }

    if (PRINT_ALL_TIME_STEPS) {
        printf("NUMR: %d\n", SIDE);
        printf("NUMC: %d\n", SIDE);
        Print_Time_Step(Temps, Now);
    }

    // Divide Temps into NUMT partitions
    struct tuple* partition_dims = Generate_Partitions();
    if (DEBUG)
        fprintf(stderr, "   partition_dims: %d %d\n", partition_dims->rows, partition_dims->cols);

    Partition_2D_Array(partition_dims->rows, partition_dims->cols);
    if (DEBUG) {
        fprintf(stderr, "   partitioning 2d array is not causing the segfault\n");
        for(int i = 0; i < NUMT; i++)
            fprintf(stderr, "       NUMT %d : [%d, %d] [%d, %d]\n", i, partitions[i].row_start, partitions[i].row_end,
                                                                    partitions[i].col_start, partitions[i].col_end);
    }
    
    // Returns the current wall clock time in seconds
    double time_init = omp_get_wtime();

    #pragma omp parallel default(none) shared(Temps,Now,Next) 
    {
        int me = omp_get_thread_num();
        DoAllWork(me);
    }

    // Calculate time elapsed from start of function until all threads completed
    double time_end = omp_get_wtime();
    double usecs = 1000000 * (time_end - time_init);
    double mega_elem_per_sec = (float)NUM_TIME_STEPS * (float)NUME / usecs;

    if (VERIFY_RESULTS) {
        printf("Performance in MegaElements/s: %10.2lf\n", mega_elem_per_sec);

        // Verify the final temperatures sum to 100
        float sum = 0;
        for (int i = 0; i < SIDE; i++) {
            for (int j = 0; j < SIDE; j++) {
                sum += Temps[Now][i][j];
            }
        }
        printf("final sum %.2f\n", sum);
    }

    if (CSV)
        printf("%2d, %8d, %10.2lf\n", NUMT, NUME, mega_elem_per_sec);
}

void DoAllWork(int me) {
    for (int step = 0; step < NUM_TIME_STEPS; step++) {

        for (int row = partitions[me].row_start; row <= partitions[me].row_end; row++) {

            // leftmost element at col_start
            {
                int col = partitions[me].col_start;

                float left = col == 0 ? 0 : Temps[Now][row][col-1];

                float right = Temps[Now][row][col+1];

                float up = row == 0? 0 : Temps[Now][row-1][col];

                float down = row == SIDE-1 ? 0 : Temps[Now][row+1][col];

                Temps[Next][row][col] = Temps[Now][row][col] + CALC_DTEMP(Temps[Now][row][col], left, right, up, down);
            }

            // middle elements: (col_start, col_end)
            for (int col = partitions[me].col_start + 1; col < partitions[me].col_end; col++) {
                
                float left = Temps[Now][row][col-1];
                float right = Temps[Now][row][col+1];

                float up = row == 0 ? 0 : Temps[Now][row-1][col];
                
                float down = row == SIDE-1 ? 0 : Temps[Now][row+1][col];
                
                Temps[Next][row][col] = Temps[Now][row][col] + CALC_DTEMP(Temps[Now][row][col], left, right, up, down);
            }

            // rightmost element at col_end
            {
                int col = partitions[me].col_end;

                float left = Temps[Now][row][col-1];

                float right = col == SIDE-1 ? 0 : Temps[Now][row][col+1];

                float up = row == 0 ? 0 : Temps[Now][row-1][col];

                float down = row == SIDE-1 ? 0 : Temps[Now][row+1][col];

                Temps[Next][row][col] = Temps[Now][row][col] + CALC_DTEMP(Temps[Now][row][col], left, right, up, down);
            }
        }

        // All threads need to wait here so that all Next values are filled
        #pragma omp barrier

        // Switch Now and Next
        #pragma omp single
        {
            Now = Next;
            Next = 1 - Next;

            if (VERIFY_RESULTS) {
                printf("\nTime step: %i\n", step);
                    
                for (int i = 0; i < SIDE; i++) {
                    for (int j = 0; j < SIDE; j++) {
                        printf(" %.2f ", Temps[Now][i][j]);
                    }
                    printf("\n");
                }
            }

            if (PRINT_ALL_TIME_STEPS) {
                Print_Time_Step(Temps, Now);
            }
        }
    }    
}