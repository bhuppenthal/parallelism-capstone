/*
    Heat Diffusion Equation Simulation - 2D Using OpenMP

    Naive implementation using vertical partitions - traversing through paritions in row-major order.

*/

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

#include "partition.h"

#ifndef PRINT_LAST_TIME_STEP
#define PRINT_LAST_TIME_STEP        false        // set to true to allow only the final time step to print
#endif

#ifndef CSV
#define CSV                         false       // set to true to print CSV of performances
#endif

float   Temps[2][SIDE][SIDE];

int     Now;                                    // which array is the "current values" = 0 or 1
int     Next;                                   // which array is being filled = 1 or 0

const int PARTITION_ROWS = 1;
const int PARTITION_COLS = 4;

void    DoAllWork(int);

int main(void) {

    Now = 0;
    Next = 1;

    // setting up the number of threads that will be used
    omp_set_num_threads(NUMT);

    for (int i = 0; i < SIDE; i++)
        for (int j = 0; j < SIDE; j++)
            Temps[Now][i][j] = (float) 0;
    
    Temps[Now][SIDE/2][SIDE/2] = (float) 100;

    if (PRINT_ALL_TIME_STEPS) {
        for (int i = 0; i < SIDE; i++) {
            for (int j = 0; j < SIDE; j++)
                printf(" %.2f ", Temps[Now][i][j]);
            printf("\n");
        }
    }

    // Calling partition function to break up sections that each thread is responsible for
    Partition_2D_Array(PARTITION_ROWS, PARTITION_COLS);

    // Returns the current wall clock time in seconds
    double time_init = omp_get_wtime();

    #pragma omp parallel default(none) shared(Temps,Now,Next) 
    {
        // save the thread number
        int me = omp_get_thread_num();
        // each thread calls this function passing in their number
        DoAllWork(me);
    }

    // Calculate time elapsed from start of function until all threads completed
    double time_end = omp_get_wtime();
    double usecs = 1000000 * (time_end - time_init);
    double mega_elem_per_sec = (float)NUM_TIME_STEPS * (float)NUME / usecs;


    if (CSV) {

        fprintf(stderr, "%2d, %8d, %10.2lf\n", NUMT, NUME, mega_elem_per_sec);

    } else {

        if (PRINT_LAST_TIME_STEP) {
            printf("Time Step: %i\n", NUM_TIME_STEPS - 1);

            for (int i = 0; i < SIDE; i++) {
                for (int j = 0; j < SIDE; j++) {
                    printf(" %.2f ", Temps[Now][i][j]);
                }
                printf("\n");
            }
        }

        fprintf(stderr, "Performance in MegaElements/s: %10.2lf\n", mega_elem_per_sec);

        // Verify the final temperatures sum to 100
        float sum = 0;
        for (int i = 0; i < SIDE; i++) {
            for (int j = 0; j < SIDE; j++) {
                sum += Temps[Now][i][j];
            }
        }

        printf("final sum %.2f\n", sum);
    }


}


void DoAllWork(int me) {

    // Determine which columns and rows this thread calculates by checking partition struct
    int first_col = partitions[me].col_start;
    int last_col = partitions[me].col_end;
    int first_row = partitions[me].row_start;
    int last_row = partitions[me].row_end;


    printf("I am thread %i\n", me);
    printf("first col is %i, last col is %i\n", first_col, last_col);
    printf("first row is %i, last row is %i\n", first_row, last_row);



    for (int step = 0; step < NUM_TIME_STEPS; step++) {

        // for each row this thread is responsible for
        for (int row = first_row; row <= last_row; row++) {
            
            // leftmost elements 
            {

                float left = 0;
                if (first_col != 0)
                    left = Temps[Now][row][first_col - 1];

                float right = 0;
                if (first_col != SIDE-1)
                    right = Temps[Now][row][first_col + 1];

                float up = 0;
                if (row != 0)
                    up = Temps[Now][row-1][first_col];

                float down = 0;
                if (row != SIDE - 1)
                    down = Temps[Now][row+1][first_col];
                
                Temps[Next][row][first_col] = Temps[Now][row][first_col] + CALC_DTEMP(Temps[Now][row][first_col], left, right, up, down);


            }

            // middle elements for each col the thread is responsible for

            for (int col = first_col + 1; col < last_col; col++) {
            // for (int col = 1; col < SIDE-1; col++)  {

                float left = Temps[Now][row][col-1];
                float right = Temps[Now][row][col+1];
                
                float up = 0;
                if (row != 0)
                    up = Temps[Now][row-1][col];

                float down = 0;
                if (row != SIDE-1)
                    down = Temps[Now][row+1][col];
                

                Temps[Next][row][col] = Temps[Now][row][col] + CALC_DTEMP(Temps[Now][row][col], left, right, up, down);

            }

            // rightmost elements
            {
                float left = 0;
                if (last_col != 0)
                    left = Temps[Now][row][last_col-1];
    
                float right = 0;
                if (last_col != SIDE-1)
                    right = Temps[Now][row][last_col+1];

                float up = 0;
                if (row != 0)
                    up = Temps[Now][row-1][last_col];

                float down = 0;
                if (row != SIDE-1)
                    down = Temps[Now][row+1][last_col];

                Temps[Next][row][last_col] = Temps[Now][row][last_col] + CALC_DTEMP(Temps[Now][row][last_col], left, right, up, down);
            }

        }


        //all threads need to wait here until all the Temps values are filled
        #pragma omp barrier

        // switching now and next and printing time steps if indicated
        #pragma omp single
        {
            Now = Next;
            Next = 1 - Next;

            if (PRINT_ALL_TIME_STEPS) {
                printf("Time Step: %i\n", step);

                for (int i = 0; i < SIDE; i++) {
                    for (int j = 0; j < SIDE; j++) {
                        printf(" %.2f ", Temps[Now][i][j]);
                    }
                    printf("\n");
                }
            }

        }
        
    }
}