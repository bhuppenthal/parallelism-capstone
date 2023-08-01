#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "heat.h"
#include "partition.h"

#define BOSS 0

#ifndef GRID_SIZE
#define GRID_SIZE 64    // global 2D array size. Total number of elements = GRID_SIZE X GRID_SIZE
#endif

#define NUMELEMENTS (GRID_SIZE * GRID_SIZE)
#define NUM_TIME_STEPS 4
#define DEBUG true
//#define WANT_EACH_TIME_STEPS_DATA

int     NumCpus; // total # of cpus involved

// number of elements in local 2D partition is PPRows x PPCols
int     PPRows;  // per-processor local partition length; the vertical partition length is equal to the GRID_SIZE
int     PPCols;  // per-processor local partition width;  the vertical partition width is equal to GRID_SIZE/NumCpus

float** PPTemps;   // per-processor local 2d array temperature data
float** NextTemps; // per-processor 2d array to hold next-values
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data

void DoOneTimeStep(int me, struct tuple* partition_dims);
// void GatherResult(int me);


int main(int argc, char *argv[]) {

    MPI_Init( &argc, &argv );
    MPI_Status status;

    int me; // which one I am - the rank of the processor

    MPI_Comm_size( MPI_COMM_WORLD, &NumCpus );
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    // Dividing the 2D array into NUMCPU partitions
    struct tuple* partition_dims = Generate_Partitions();
    
    // DEBUGGING
    // printf("My name is %d: I have %d rows and %d cols\n", me, partition_dims->rows, partition_dims->cols);
    // DEBUGGING

    // what is the length (PPCols) and width (PPRow) of the processor's rectangular partition
    Partition_2D_Array(partition_dims->rows, partition_dims->cols);
    PPRows = partitions[me].row_end - partitions[me].row_start + 1;
    PPCols = partitions[me].col_end - partitions[me].col_start + 1;

    // DEBUGGING
    // printf("My name is rank %i: My row start is %d and my row end is %d\n", me, partitions[me].row_start, partitions[me].row_end);
    // printf("My name is rank %i: My col start is %d and my col end is %d\n", me, partitions[me].col_start, partitions[me].col_end);
    //DEBUGGING

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

        //DEBUGGING - PRINT BOSS PP TEMPS
        // printf("Hi I am BOOS my PPRows is %d and my PPCols is %d\n", PPRows, PPCols);
        // printf("Hi I am BOSS here are my PP Temps:\n");
        // for (int i = 0; i < PPRows; i++) {
        //     for (int j = 0; j < PPCols; j++) {
        //         printf("%2f  ",PPTemps[i][j]);
        //     }
        //     printf("\n");
        // }
        // DEBUGGING

        // BOSS will send the rest of the partitions to the other processors sending a “row segment” 
        for (int dest = 1; dest < NumCpus; dest++) {
            int startRow = partitions[dest].row_start;
            int endRow = partitions[dest].row_end; 
            int rowIdx = startRow; // Index of the next row to be sent
            int col_start = partitions[dest].col_start;

            while (rowIdx <= endRow) {
                MPI_Send(&TempData[rowIdx][col_start], PPCols, MPI_FLOAT, dest, rowIdx, MPI_COMM_WORLD);
                rowIdx++;
            }
        }
    } else {
        // Receive the plate data from the BOSS Processor
        int rowIdx = partitions[me].row_start; //Index of next row to be received

        for (int i = 0; i < PPRows; i++) {
            MPI_Recv(&PPTemps[i][0], PPCols, MPI_FLOAT, BOSS, rowIdx, MPI_COMM_WORLD, &status);
            rowIdx++;
        }

        //DEBUGGING
        // printf("Hi I am rank %i my PPRows is %d and my PPCols is %d\n", me, PPRows, PPCols);
        // printf("Hello my rank is %i, and here is my PPTemps:\n", me);
        // for (int i = 0; i < PPRows; i++) {
        //     for (int j = 0; j < PPCols; j++) {
        //         printf("%2f  ",PPTemps[i][j]);
        //     }
        //     printf("\n");
        // }
        //DEBUGGING

    }

    // all the vertical paritions (ie. PPTemps 2D arrays) have been filled for each processor and BOSS
    double time0 = MPI_Wtime();

    for(int steps = 0; steps < NUM_TIME_STEPS; steps++) {
        // do the computation for one time step:
        DoOneTimeStep(me, partition_dims);
        // ask for all the data:
#ifdef WANT_EACH_TIME_STEPS_DATA
        // GatherResult(me);
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
    // GatherResult(me);
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


void DoOneTimeStep(int me, struct tuple* partition_dims) {
    
    MPI_Status status;

    // 1D arrays that contain the border left, right, up, and bottom boundary elements of each 
    // should I change these to arrays?? 
    // PPTempsLeft = new float[PPRows];
    // PPTempsRight = new float[PPRows];
    // PPTempsTop = new float[PPCols];
    // PPTempsBottom = new float[PPCols];

    std::vector<float> PPTempsLeft(PPRows, 0.0);
    std::vector<float> PPTempsRight(PPRows, 0.0);
    std::vector<float> PPTempsTop(PPCols, 0.0);
    std::vector<float> PPTempsBottom(PPCols, 0.0);

    
    // send out left boundary values
    // (tag is point of view of the sender)
     if (partitions[me].col_start != 0) // if i am not one of the leftmost groups
    {
        // Need to copy all left boundary elements into PPTempsLeft array to send them
        for (int i = 0; i < PPRows; i++) {
            PPTempsLeft[i] = PPTemps[i][0];
        }

        // send PPTempsLeft[0] temps to me-1 using the ‘L’ tag
        MPI_Send(&PPTempsLeft[0], PPRows, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD);
        if(DEBUG) {
            fprintf(stderr, "%3d sent 'L' to %3d\n", me, me-1);
            // printf("PPTemps Left for rank %i: \n", me);
            // for (int i =0; i < PPRows; i++) {
            //     printf("%0.2f\n", PPTempsLeft[i]);
            // }
        }
    }  

    // send out the right boundary values 
    // (tag is point of view of the sender)
    if (partitions[me].col_end != GRID_SIZE-1) //If I am not one of the rightmost groups
    {
        // Need to copy all right boundary elements into PPTempsRight array to send them
	    for (int i = 0; i < PPRows; i++) {
	        PPTempsRight[i] = PPTemps[i][PPCols-1];
        }

        //send PPTempsRight[0] temps - to me+1 using the ‘R’ tag 
	    MPI_Send(&PPTempsRight[0], PPRows, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD);
	    if(DEBUG) fprintf(stderr, "%3d sent 'R' to %3d\n", me, me+1);
    }

    // send out top row boundary elements 
    if (partitions[me].row_start != 0) // if I'm not one of the groups on the top
    {                     
        // send my PPTemps[0] to partition below useing tag 'T'
        // me - partition_dims->cols
        MPI_Send(PPTemps[0], PPCols, MPI_FLOAT, me - partition_dims->cols, 'T', MPI_COMM_WORLD);
        if(DEBUG) fprintf(stderr, "%3d sent 'T' to %3d\n", me, me - partition_dims->cols);
    }

    // send out bottom row boundary elements
    if (partitions[me].row_end != GRID_SIZE-1) // if i'm not one of the groups on the bottom
    {
        // send my PPTemps[PPRows - 1] to partition above using tag 'B'
        // me + partition_dims->cols?
        MPI_Send(PPTemps[PPRows - 1], PPCols, MPI_FLOAT, me + partition_dims->cols, 'B', MPI_COMM_WORLD);
        if(DEBUG) fprintf(stderr, "%3d sent 'B' to %3d\n", me, me + partition_dims->cols);
    }

    // recieve all the boundary arrays
    std::vector<float> left(PPRows, 0.0);
    std::vector<float> right(PPRows, 0.0);
    std::vector<float> up(PPCols, 0.0);
    std::vector<float> down(PPCols, 0.0);

    if (partitions[me].col_start != 0)   // if i'm not the first partition on the left
    {              
        // receive my "left" from me - 1 using tag 'R'
        MPI_Recv(&left[0], PPRows, MPI_FLOAT, me - 1, 'R', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf( stderr, "%3d received 'R' from %3d\n", me, me - 1);
    }

    if (partitions[me].col_end != GRID_SIZE -1) // if not the last partition on the right
    {        
        // receive my "right" from me+1 using tag 'L'
        MPI_Recv(&right[0], PPRows, MPI_FLOAT, me + 1, 'L', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf(stderr, "%3d received 'L' from %3d\n", me, me + 1);
    }

    if (partitions[me].row_start != 0) // if I'm not one of the groups on the top
    {
        // recieve my "up" from from me - partition_dim->cols using tag "B"
        MPI_Recv(&up[0], PPCols, MPI_FLOAT, me - partition_dims->cols, 'B', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf(stderr, "%3d received 'B' from %3d\n", me, me - partition_dims->cols);
    }

    if (partitions[me].row_end != GRID_SIZE-1) // if i'm not one of the groups on the bottom
    {
        // recieve my "down" from me + parititon-dims->cols using tag "T"
        MPI_Recv(&down[0], PPCols, MPI_FLOAT, me + partition_dims->cols, 'T', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf(stderr, "%3d received 'T' from %3d\n", me, me + partition_dims->cols);
    }

}