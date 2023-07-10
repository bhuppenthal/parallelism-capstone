/*
    Heat Diffusion Equation Simulation - 2D Using OpenMP

    Naive implementation using horizontal partitions.
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
#define PRINT_ALL_TIME_STEPS		true         // set to true to allow all time steps to print
#endif                            

// TODO: add a CACHE_ALIGNMENT constant
// 0 to 63, experiment with various misalignments to see results
// with single precision floats, 16 floats / cache line

const int SIDE = (int) sqrt(NUME);

// TODO: switch back over to a flip flopping now / next buffer (1D example)
float   Temps[2][SIDE][SIDE];

int     Now = 0;                                    // which array is the "current values" = 0 or 1
int     Next = 1;                                   // which array is being filled = 1 or 0

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

    // Setting number of threads that will be used
    omp_set_num_threads(NUMT);

    // Setting all initial temperatures to 0, except the middle value which is set to 100
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

    float sum = 0;
    for (int i = 0; i < SIDE; i++) {
        for (int j = 0; j < SIDE; j++) {
            sum += Temps[Now][i][j];
        }
    }

    printf("final sum %.2f\n", sum);
}

void DoAllWork(int me) {
    // Determine what rows this thread calculates
    int first_row = me * SIDE / NUMT;
    int last_row = first_row + (SIDE / NUMT - 1);

    for (int step = 0; step < NUM_TIME_STEPS; step++) {

        // for each row assigned to this thread
        for (int i = first_row; i <= last_row; i++) {

            // leftmost element
            {
                float left = 0;
                float right = Temps[Now][i][1];

                float up = 0;
                if (i != 0)
                    up = Temps[Now][i-1][0];

                float down = 0;
                if (i != SIDE - 1)
                    down = Temps[Now][i+1][0];

                Temps[Next][i][0] = Temps[Now][i][0] + CALC_DTEMP(Temps[Now][i][0], left, right, up, down);
            }

            // middle elements
            for (int j = 1; j < SIDE - 1; j++) {
                float up = 0;
                if (i != 0)
                    up = Temps[Now][i-1][j];
                
                float down = 0;
                if (i != SIDE - 1)
                    down = Temps[Now][i+1][j];
                
                float left = Temps[Now][i][j-1];
                float right = Temps[Now][i][j+1];

                Temps[Next][i][j] = Temps[Now][i][j] + CALC_DTEMP(Temps[Now][i][j], left, right, up, down);
            }

            // rightmost element
            {
                float left = Temps[Now][i][SIDE-2];
                float right = 0;

                float up = 0;
                if (i != 0)
                    up = Temps[Now][i-1][SIDE-1];

                float down = 0;
                if (i != SIDE - 1)
                    down = Temps[Now][i+1][SIDE-1];

                Temps[Next][i][SIDE-1] = Temps[Now][i][SIDE-1] + CALC_DTEMP(Temps[Now][i][SIDE-1], left, right, up, down);
            }
            
        }

        // all threads need to wait here so that all NowTemps values are filled:
        #pragma omp barrier

        // Switch Now and Next.
        Now = Next;
        Next = 1 - Next;

        if (PRINT_ALL_TIME_STEPS) {
            #pragma omp single
            {
                printf("\nTime step: %i\n", step);
                
                for (int i = 0; i < SIDE; i++) {
                    for (int j = 0; j < SIDE; j++) {
                        printf(" %.2f ", Temps[Now][i][j]);
                    }
                    printf("\n");
                }
            } // implied barrier exists here
        }
    }    
}