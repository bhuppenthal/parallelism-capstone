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

int     now;                                    // which array is the "current values" = 0 or 1
int     next;                                   // which array is being filled = 1 or 0

// Heat Diffusion Equation Constants
const float RHO = 8050;
const float C = 0.466;
const float K = 20;
// float k_over_rho_c = K/(RHO*C);                 // units of m^2/s NOTE: Cannot be a const (true for OpenMP too?)
// // K/(RHO*C) = 5.33 x 10^-6 m^2/s

const float DX = 1;
const float DT = 1;

void    DoAllWork(int);

int main(void) {

    // Setting number of threads that will be used
    omp_set_num_threads(NUMT);

    now = 0;
    next = 1;

    // Setting all initial temperatures to 0, except the middle value which is set to 100
    for (int i = 0; i < NUME; i++)
        Temps[now][i] = 0;
    
    Temps[now][NUME/2] = 100;

    // Returns the current wall clock time in seconds
    double time_init = omp_get_wtime();

    #pragma omp parallel default(none) shared(Temps,now,next) 
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

    
}

void DoAllWork(int me) {

    // What range of global Temps array this thread is responsible for:
    int first = me * NUM_ELEM_PER_THREAD;
    int last = first + (NUM_ELEM_PER_THREAD - 1);

    for (int step = 0; step < NUM_TIME_STEPS; step++) {
        
        // first element on the left:
        {
            float left = 0;
            if (me != 0)
                left = Temps[now][first - 1];
            
            float dtemp = ((K/(RHO*C)) * 
                        (left - (2.0 * Temps[now][first]) + Temps[now][first + 1]) / (DX * DX)) * DT;
            
            Temps[next][first] = Temps[now][first] + dtemp;
            

        }
    }    

}