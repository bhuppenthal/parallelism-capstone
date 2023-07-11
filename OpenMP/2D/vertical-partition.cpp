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
#define PRINT_ALL_TIME_STEPS		false       // set to true to allow all time steps to print
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

}


void DoAllWork(int me) {
    
}