/*

    Heat Diffusion Equation Simulation - 1D Using OpenMP
    Code Obtained From:
    Bailey, Mike. (2023). Data Decomposition Slides 6 - 9. Retrieved on 07/05/2023 From:
    https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/data_decomposition.1pp.pdf 

*/

#include <stdio.h>
#include <math.h>
#include <omp.h>

#define NUM_TIME_STEPS  100                     // number of time steps the simultation runs through

#ifndef NUME
#define NUME            1024                    // total number of elements
#endif 

#ifndef NUMT
#define NUMT            4                       // number of threads to use
#endif

#define NUM_ELEM_PER_THREAD    (NUME/NUMT)      // number of elements in each thread

float   Temps[2][NUME];                         // storing all temperatures "Now" and "Next" states

int     Now;                                    // which array is the "current values" = 0 or 1
int     Next;                                   // which array is being filled = 1 or 0

void    DoAllWork(int);

int main(void) {

    printf("Hello world\n");
    
}

void DoAllWork(int me) {
    

}