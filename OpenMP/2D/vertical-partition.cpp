/*
    Heat Diffusion Equation Simulation - 2D Using OpenMP

    Naive implementation using vertical partitions.

*/

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

#define NUM_TIME_STEPS  100                      // number of time steps the simultation runs through

#ifndef NUME
#define NUME            64                    // total number of elements
#endif

#ifndef NUMT
#define NUMT            4                       // number of threads to use
#endif

#define NUM_ELEM_PER_THREAD    (NUME/NUMT)      // number of elements in each thread

#ifndef PRINT_ALL_TIME_STEPS
#define PRINT_ALL_TIME_STEPS		true       // set to true to allow all time steps to print
#endif

#ifndef PRINT_LAST_TIME_STEP
#define PRINT_LAST_TIME_STEP        false        // set to true to allow only the final time step to print
#endif

#ifndef CSV
#define CSV                         false       // set to true to print CSV of performances
#endif                          

const int SIDE = (int) sqrt(NUME);

float   Temps[2][SIDE][SIDE];

int     Now;                                    // which array is the "current values" = 0 or 1
int     Next;                                   // which array is being filled = 1 or 0

// Heat Diffusion Equation Constants
const float RHO = 8050;
const float C = 0.466;
const float K = 20;
float k_over_rho_c = K/(RHO*C);                 // units of m^2/s NOTE: Cannot be a const (true for OpenMP too?)
// K/(RHO*C) = 5.33 x 10^-6 m^2/s

const float DX = 1;
const float DY = 1;
const float DT = 1;

// Heat Diffusion Equation
#define CALC_DTEMP(elem, left, right, up, down) (k_over_rho_c * DT * ((down - 2*elem + up)/(DY * DY) + (left - 2*elem + right)/(DX * DX)));

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


void DoAllWork(int me) {

    // Determine which columns this thread calculates
    int first_col = (SIDE / NUMT) * me;
    int last_col = first_col + ((SIDE / NUMT) - 1);

    // printf("First col: %i, last col: %i\n", first_col, last_col);

    for (int step = 0; step < NUM_TIME_STEPS; step++) {

        // printf("STEP: %d\n", step);

        // for each column this thread is responsible for
        for (int col = first_col; col <= last_col; col++) {
            
            // printf("I am thread %d, this is col: %d\n", me, col);

            // uppermost elements 
            {
                float up = 0;
                float down = Temps[Now][1][col];

                float left = 0;
                if (col != 0) 
                    left = Temps[Now][0][col-1];
                
                float right = 0;
                if (col != SIDE -1) 
                    right = Temps[Now][0][col+1];
                
                Temps[Next][0][col] = Temps[Now][0][col] + CALC_DTEMP(Temps[Now][0][col], left, right, up, down);

            }


            // middle elements
            for (int row = 1; row <= SIDE - 2; row++) {
                
                float left = 0;
                if (col != 0) 
                    left = Temps[Now][row][col-1];
                

                float up = Temps[Now][row-1][col];
                float down = Temps[Now][row+1][col];

                float right = 0;
                if (col != SIDE-1) 
                    right = Temps[Now][row][col+1];
                

                Temps[Next][row][col] = Temps[Now][row][col] + CALC_DTEMP(Temps[Now][row][col], left, right, up, down);

            }

            // bottomost elements
            {
                float up = Temps[Now][SIDE-2][col];
                float down = 0;

                float left = 0;
                if (col != 0)
                    left = Temps[Now][SIDE-1][col-1];
                
                float right = 0;
                if (col != SIDE-1)
                    right = Temps[Now][SIDE-1][col+1];

                Temps[Next][SIDE-1][col] = Temps[Now][SIDE-1][col] + CALC_DTEMP(Temps[Now][SIDE-1][col], left, right, up, down);
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