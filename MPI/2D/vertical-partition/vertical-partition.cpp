#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "heat.h"

#define BOSS 0

#ifndef GRID_SIZE
#define GRID_SIZE 8     // global 2D array size. Total number of elements = GRID_SIZE X GRID_SIZE
#endif

#define NUMELEMENTS (GRID_SIZE * GRID_SIZE)
#define NUM_TIME_STEPS 4
#define DEBUG false
//#define WANT_EACH_TIME_STEPS_DATA

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
    printf("hello world vertical partition");
}