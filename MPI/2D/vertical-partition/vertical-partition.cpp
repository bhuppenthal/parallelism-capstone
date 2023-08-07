#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "heat.h"

#define BOSS 0

#ifndef GRID_SIZE
#define GRID_SIZE 32    // global 2D array size. Total number of elements = GRID_SIZE X GRID_SIZE
#endif

#define NUMELEMENTS (GRID_SIZE * GRID_SIZE)
#define NUM_TIME_STEPS 100
#define DEBUG false
// #define WANT_EACH_TIME_STEPS_DATA

int     NumCpus; // total # of cpus involved

// number of elements in local 2D partition is PPRows x PPCols
int     PPRows;  // per-processor local partition length; the vertical partition length is equal to the GRID_SIZE
int     PPCols;  // per-processor local partition width;  the vertical partition width is equal to GRID_SIZE/NumCpus

float** PPTemps;   // per-processor local 2d array temperature data
float** NextTemps; // per-processor 2d array to hold next-values
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data

void DoOneTimeStep(int);
void GatherResult(int me);

int main(int argc, char *argv[]) {

    MPI_Init( &argc, &argv );
    MPI_Status status;

    int me; // which one I am - the rank of the processor

    MPI_Comm_size( MPI_COMM_WORLD, &NumCpus );
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    PPRows = GRID_SIZE;
    PPCols = GRID_SIZE / NumCpus;

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
    
        // Copy the BOSS vertical parition to BOSS local 2D array (not MPI_Send) since BOSS does will do the computations for the first partition
        for (int i = 0; i < PPRows; i++) {
            for (int j = 0; j < PPCols; j++) {
                PPTemps[i][j] = TempData[i][j];
            }
        }

        // BOSS will send the rest of the partitions to the other processors sending a “row segment” 
        for (int dest = 1; dest < NumCpus; dest++) {
            int startRow = 0;
            int endRow = PPRows - 1; 
            int rowIdx = startRow; // Index of the next row to be sent
            int col_start = PPCols * dest;

            while (rowIdx <= endRow) {
                MPI_Send(&TempData[rowIdx][col_start], PPCols, MPI_FLOAT, dest, rowIdx, MPI_COMM_WORLD);
                rowIdx++;
            }
        }
    } else {
        // Receive the plate data from the BOSS Processor
        int rowIdx = 0; //Index of next row to be received

        for (int i = 0; i < PPRows; i++) {
            MPI_Recv(&PPTemps[i][0], PPCols, MPI_FLOAT, BOSS, rowIdx, MPI_COMM_WORLD, &status);
            rowIdx++;
        }

    }

    // all the vertical paritions (ie. PPTemps 2D arrays) have been filled for each processor and BOSS
    double time0 = MPI_Wtime();

    for(int steps = 0; steps < NUM_TIME_STEPS; steps++) {
        // do the computation for one time step:
        DoOneTimeStep(me);
        // ask for all the data:
#ifdef WANT_EACH_TIME_STEPS_DATA
        GatherResult(me);
        //MPI_Barrier(MPI_COMM_WORLD);

        if (me == BOSS) {
            fprintf(stdout, "Time step: %3d\n", steps);
            if (me == BOSS) {
                for (int i = 0; i < GRID_SIZE; i++) {
                    for (int j = 0; j < GRID_SIZE; j++) {
                        std::cout << std::fixed << std::setprecision(2) << TempData[i][j] << " ";
                    }
                    std::cout << std::endl;
                }            
            }
        }
#endif    
    }
#ifndef WANT_EACH_TIME_STEPS_DATA
    GatherResult(me);
#endif

    double time1 = MPI_Wtime( );

    if( me == BOSS )
    {
        double seconds = time1 - time0;
        double performance = 
            (double)NUM_TIME_STEPS * (double)NUMELEMENTS / seconds / 1000000.;
                // mega-elements computed per second
        fprintf( stderr, "%3d, %10d, %8.2lf\n", NumCpus, NUMELEMENTS, performance );
    }

    // Deallocate memory for PPTemps
    for (int i = 0; i < PPRows; i++) {
        delete[] PPTemps[i];
    }
    delete[] PPTemps;

    // Deallocate memory for NextTemps
    for (int i = 0; i < PPRows; i++) {
        delete[] NextTemps[i];
    }
    delete[] NextTemps;

    // Deallocate memory for TempData (only if BOSS)
    if (me == BOSS) {
        for (int i = 0; i < GRID_SIZE; i++) {
            delete[] TempData[i];
        }
        delete[] TempData;
    }

    MPI_Finalize( );
    return 0;

}

void DoOneTimeStep(int me) {
    
    MPI_Status status;

    float PPTempsLeft[GRID_SIZE] = {0.};
    float PPTempsRight[GRID_SIZE] = {0.};
    
    // send out left and right end values
    // (tag is point of view of the sender)

    if (me != 0) // if i am not the first group on the left
    {
        // Need to copy all left boundary elements into PPTempsLeft array to send them
        for (int i = 0; i < PPRows; i++) {
            PPTempsLeft[i] = PPTemps[i][0];
        }
        // send PPTempsLeft[0] temps to me-1 using the ‘L’ tag
        MPI_Send(&PPTempsLeft[0], PPRows, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD);
        if(DEBUG) {
            fprintf(stderr, "%3d sent 'L' to %3d\n", me, me-1);
        }
    }

    if (me != NumCpus - 1) // if I am not the last group on the right
    {
        // Need to copy all right boundary elements into PPTempsRight array to send them
	    for (int i = 0; i < PPRows; i++) {
	        PPTempsRight[i] = PPTemps[i][PPCols-1];
        }

        //send PPTempsRight[0] temps - to me+1 using the ‘R’ tag 
	    MPI_Send(&PPTempsRight[0], PPRows, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD);
	    if(DEBUG) fprintf(stderr, "%3d sent 'R' to %3d\n", me, me+1);
    }

    // Recieve left and right boundary values into left and right arrays respectively
    float left[GRID_SIZE] = {0.};
    float right[GRID_SIZE] = {0.};

    if (me != 0)   // if i'm not the first partition on the left
    {              
        // receive my "left" from me - 1 using tag 'R'
        MPI_Recv(&left[0], GRID_SIZE, MPI_FLOAT, me - 1, 'R', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf( stderr, "%3d received 'R' from %3d\n", me, me - 1);
    }

    if (me != NumCpus - 1) // if not the last partition on the right
    {        
        // receive my "right" from me+1 using tag 'L'
        MPI_Recv(&right[0], GRID_SIZE, MPI_FLOAT, me + 1, 'L', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf(stderr, "%3d received 'L' from %3d\n", me, me + 1);
    }

    // Now all processors have the necessary boundary elements we can calculate the NextTemps of each vertical partition
    float up = 0.;
    float down = 0.;

    // first row on the top
    // top-left corner element
    NextTemps[0][0] = PPTemps[0][0] + CALC_DTEMP(PPTemps[0][0], left[0], PPTemps[0][1], up, PPTemps[1][0]);


    // middle elements in the first row
    for (int i = 1; i < PPCols - 1; i++) {
        NextTemps[0][i] = PPTemps[0][i] + CALC_DTEMP(PPTemps[0][i], PPTemps[0][i - 1], PPTemps[0][i + 1], up, PPTemps[1][i]);
    }


    // top-right corner element
    NextTemps[0][PPCols - 1] = PPTemps[0][PPCols - 1] + CALC_DTEMP(PPTemps[0][PPCols - 1], PPTemps[0][PPCols - 2], right[0], up, PPTemps[1][PPCols - 1]);


    // all the rows in the middle
    for (int i = 1; i < PPRows - 1; i++) {
        // left-most elements
        NextTemps[i][0] = PPTemps[i][0] + CALC_DTEMP(PPTemps[i][0], left[i], PPTemps[i][1], PPTemps[i - 1][0], PPTemps[i + 1][0]);


        for (int j = 1; j < PPCols - 1; j++) {
            NextTemps[i][j] = PPTemps[i][j] + CALC_DTEMP(PPTemps[i][j], PPTemps[i][j - 1], PPTemps[i][j + 1], PPTemps[i - 1][j], PPTemps[i + 1][j]);
        }


        // rightmost elements
        NextTemps[i][PPCols - 1] = PPTemps[i][PPCols - 1] + CALC_DTEMP(PPTemps[i][PPCols - 1], PPTemps[i][PPCols - 2], right[i], PPTemps[i - 1][PPCols - 1], PPTemps[i + 1][PPCols - 1]);
    }


    // last row on the bottom
    // bottom-left corner element
    NextTemps[PPRows - 1][0] = PPTemps[PPRows - 1][0] + CALC_DTEMP(PPTemps[PPRows - 1][0], left[PPRows - 1], PPTemps[PPRows - 1][1], PPTemps[PPRows - 2][0], down);


    // middle elements in the last row
    for (int i = 1; i < PPCols - 1; i++) {
        NextTemps[PPRows - 1][i] = PPTemps[PPRows - 1][i] + CALC_DTEMP(PPTemps[PPRows - 1][i], PPTemps[PPRows - 1][i - 1], PPTemps[PPRows - 1][i + 1], PPTemps[PPRows - 2][i], down);
    }


    // bottom-right corner element
    NextTemps[PPRows - 1][PPCols - 1] = PPTemps[PPRows - 1][PPCols - 1] + CALC_DTEMP(PPTemps[PPRows - 1][PPCols - 1], PPTemps[PPRows - 1][PPCols - 2], right[PPRows-1], PPTemps[PPRows - 2][PPCols - 1], down);


    // update the local dataset
    for (int i = 0; i < PPRows; i++) {
        for (int j = 0; j < PPCols; j++) {
            PPTemps[i][j] = NextTemps[i][j];
        }
    }

}

void GatherResult(int me) {
    MPI_Status status;

    if (me == BOSS) {
        // Copy the BOSS local temp array data to global temp array
        for (int i = 0; i < PPRows; i++) {
            memcpy(TempData[i], PPTemps[i], PPCols * sizeof(float));
        }

        // Receive data from other processors
        for (int src = 1; src < NumCpus; src++) {
            int startRow = 0;// index of first row of the strip
            int endRow = PPRows - 1;
            int start_col = src * PPCols;
            int rowIdx = startRow; // Index of the next row to be received


            while (rowIdx <= endRow) {
                MPI_Recv(&TempData[rowIdx][start_col], PPCols, MPI_FLOAT, src, rowIdx, MPI_COMM_WORLD, &status);
                if(DEBUG) fprintf(stderr, "BOSS received %3d from %3d\n", rowIdx, src);
                rowIdx++;
            }
        }


    } else {
        // Send the local temp array PPTemps data back to the BOSS processor row by row
        int startRow = 0; // Index of the first row of the vertical partition

        for (int i = 0; i < PPRows; i++) {
            int tag = startRow + i; // message tag of MPI_Send
            MPI_Send(PPTemps[i], PPCols, MPI_FLOAT, BOSS, tag, MPI_COMM_WORLD);
            if(DEBUG) fprintf(stderr, "%3d sent %3d to BOSS\n", me, tag);
        }

    }
}