#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
#include <iomanip>
// #include "heat.h"

#define BOSS 0

#ifndef GRID_SIZE
#define GRID_SIZE 8     // global 2D array size. Total number of elements = GRID_SIZE X GRID_SIZE
#endif

#define NUMELEMENTS (GRID_SIZE * GRID_SIZE)
#define NUM_TIME_STEPS 4
#define DEBUG false
// #define WANT_EACH_TIME_STEPS_DATA

int     NumCpus; // total # of cpus involved

// number of elements in local 2D partition is PPRows x PPCols
int     PPRows;  // per-processor local partition length; the vertical partition length is equal to the GRID_SIZE
int     PPCols;  // per-processor local partition width;  the vertical partition width is equal to GRID_SIZE/NumCpus

float** PPTemps;   // per-processor local 2d array temperature data
float** NextTemps; // per-processor 2d array to hold next-values
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data

float* PPTempsLeft; // per-processor local 1D array holding all values in the leftmost column for sending to other processors
float* PPTempsRight; // per-processor local 1D array holding all values in the rightmost column for sending to other processors

void DoOneTimeStep(int);
void GatherResult(int me);

int main(int argc, char *argv[]) {

    MPI_Init( &argc, &argv );

    int me; // which one I am - the rank of the processor

    MPI_Comm_size( MPI_COMM_WORLD, &NumCpus );
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    PPRows = GRID_SIZE;
    PPCols = GRID_SIZE / NumCpus;

    // the arrays to hold boundary elements that are sent between processors
    PPTempsLeft = new float[PPRows];
    PPTempsRight = new float[PPRows];

    // local 2D array (ie. vertical partition) for each CPU to hold thier section of temperatures
    PPTemps = new float* [PPRows];
    for (int i = 0; i < PPRows; i++) {
        PPTemps[i] = new float[PPCols];
    }

    NextTemps = new float* [PPRows];
    for (int i = 0; i < PPRows; i++) {
        NextTemps[i] = new float[PPCols];
    }

    if (me == BOSS) {
    // Initialize and populate the 2D plate array
    TempData = new float* [GRID_SIZE];
    for (int i = 0; i < GRID_SIZE; i++) {
        TempData[i] = new float[GRID_SIZE];
    }

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            TempData[i][j] = 0.;
        }
    }
    TempData[GRID_SIZE/2][GRID_SIZE/2] = 100.;

#ifdef WANT_EACH_TIME_STEPS_DATA
        // Print out TempData to stderr
        fprintf(stdout, "Initial Matrix\n");
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                std::cout << std::fixed << std::setprecision(2) << TempData[i][j] << " ";
            }
            std::cout << std::endl;
        }
#endif  


    MPI_Finalize( );
    return 0;


    }
}