#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include "heat.h"

#define BOSS 0

#ifndef GRID_SIZE
#define GRID_SIZE 8 // global 2d array size. Total number of elements = GRID_SIZE X GRID_SIZE
#endif

#define NUMELEMENTS (GRID_SIZE * GRID_SIZE)
#define NUM_TIME_STEPS 4
#define DEBUG false
//#define WANT_EACH_TIME_STEPS_DATA

int     NumCpus; // total # of cpus involved
int     PPRows;  // per-processor local 2d array size, local number of elements = PPRows X GRID_SIZE

float** PPTemps;   // per-processor local 2d array temperature data
float** NextTemps; // per-processor 2d array to hold next-values
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data

void DoOneTimeStep(int);
void GatherResult(int me);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Status status;

    int me;      // which one I am - the rank of a processor

    MPI_Comm_size(MPI_COMM_WORLD, &NumCpus);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    PPRows = GRID_SIZE / NumCpus;

    PPTemps = new float* [PPRows];
    for (int i = 0; i < PPRows; i++) {
        PPTemps[i] = new float[GRID_SIZE];
    }

    NextTemps = new float* [PPRows];
    for (int i = 0; i < PPRows; i++) {
        NextTemps[i] = new float[GRID_SIZE];
    }

    // broadcast the constant:
    MPI_Bcast( (void *)&k_over_rho_c, 1, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

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

        // Transfer the BOSS strip data to BOSS local array (not MPI_Send)
        for (int i = 0; i < PPRows; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                PPTemps[i][j] = TempData[i][j];
            }
        }

        // Distribute the plate data to other processors row by row
        for (int dest = 1; dest < NumCpus; dest++) {
            int startRow = dest * PPRows;
            int endRow = startRow + PPRows - 1;
            int rowIdx = startRow; // Index of the next row to be sent

            while (rowIdx <= endRow) {
                MPI_Send(&TempData[rowIdx][0], GRID_SIZE, MPI_FLOAT, dest, rowIdx, MPI_COMM_WORLD);
                if(DEBUG) fprintf(stderr, "BOSS sent %3d to %3d\n", rowIdx, dest);
                rowIdx++;
            }
        }
    } else {
        // Receive the plate data from the BOSS processor
        int rowIdx = me * PPRows; // Index of the next row to be received

        for (int i = 0; i < PPRows; i++) {
            MPI_Recv(&PPTemps[i][0], GRID_SIZE, MPI_FLOAT, BOSS, rowIdx, MPI_COMM_WORLD, &status);
            if(DEBUG) fprintf( stderr, "%3d received %3d from BOSS\n", me, rowIdx);
            rowIdx++;
        }
    }

    // all the PPTemps arrays have now been filled
    // do the time steps:
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

    MPI_Finalize();
    return 0;
}

void DoOneTimeStep(int me) {
    MPI_Status status;

    // send out the top and bottom rows values
    // the tag is from the point view of the sender
    if (me != 0) {                     // if not the first strip on the top
        // send my PPTemps[0] to me-1 useing tag 'T'
        MPI_Send(PPTemps[0], GRID_SIZE, MPI_FLOAT, me - 1, 'T', MPI_COMM_WORLD);
        if(DEBUG) fprintf(stderr, "%3d sent 'T' to %3d\n", me, me - 1);
    }

    if (me != NumCpus - 1) {            // if not the last strip on the bottom
        // send my PPTemps[PPRows - 1] to me+1 using tag 'B'
        MPI_Send(PPTemps[PPRows - 1], GRID_SIZE, MPI_FLOAT, me + 1, 'B', MPI_COMM_WORLD);
        if(DEBUG) fprintf(stderr, "%3d sent 'B' to %3d\n", me, me + 1);
    }

    float up[GRID_SIZE] = {0.};
    float down[GRID_SIZE] = {0.};

    if (me != 0) {                      // if not the first strip on the top
        // receive my "up" from me-1 using tag 'B'
        MPI_Recv(&up[0], GRID_SIZE, MPI_FLOAT, me - 1, 'B', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf( stderr, "%3d received 'B' from %3d\n", me, me - 1);
    }

    if (me != NumCpus - 1) {            // if not the last strip on the bottom
        // receive my "down" from me+1 using tag 'T'
        MPI_Recv(&down[0], GRID_SIZE, MPI_FLOAT, me + 1, 'T', MPI_COMM_WORLD, &status);
        if(DEBUG) fprintf(stderr, "%3d received 'T' from %3d\n", me, me + 1);
    }


    float left = 0.;
    float right = 0.;

    // first row on the top
    // top-left corner element
    NextTemps[0][0] = PPTemps[0][0] + CALC_DTEMP(PPTemps[0][0], left, PPTemps[0][1], up[0], PPTemps[1][0]);

    // middle elements in the first row
    for (int i = 1; i < GRID_SIZE - 1; i++) {
        NextTemps[0][i] = PPTemps[0][i] + CALC_DTEMP(PPTemps[0][i], PPTemps[0][i - 1], PPTemps[0][i + 1], up[i], PPTemps[1][i]);
    }

    // top-right corner element
    NextTemps[0][GRID_SIZE - 1] = PPTemps[0][GRID_SIZE - 1] + CALC_DTEMP(PPTemps[0][GRID_SIZE - 1], PPTemps[0][GRID_SIZE - 2], right, up[GRID_SIZE - 1], PPTemps[1][GRID_SIZE - 1]);

    // all the rows in the middle
    for (int i = 1; i < PPRows - 1; i++) {
        // left-most elements
        NextTemps[i][0] = PPTemps[i][0] + CALC_DTEMP(PPTemps[i][0], left, PPTemps[i][1], PPTemps[i - 1][0], PPTemps[i + 1][0]);

        for (int j = 1; j < GRID_SIZE - 1; j++) {
            NextTemps[i][j] = PPTemps[i][j] + CALC_DTEMP(PPTemps[i][j], PPTemps[i][j - 1], PPTemps[i][j + 1], PPTemps[i - 1][j], PPTemps[i + 1][j]);

            //int temp = k_over_rho_c * DT * ((PPTemps[i + 1][j] - 2*PPTemps[i][j] + PPTemps[i - 1][j])/(DY * DY) + (PPTemps[i][j - 1] - 2*PPTemps[i][j] + PPTemps[i][j + 1])/(DX * DX));
            //NextTemps[i][j] = PPTemps[i][j] + temp;
        }

        // right-most elements
        NextTemps[i][GRID_SIZE - 1] = PPTemps[i][GRID_SIZE - 1] + CALC_DTEMP(PPTemps[i][GRID_SIZE - 1], PPTemps[i][GRID_SIZE - 2], right, PPTemps[i - 1][GRID_SIZE - 1], PPTemps[i + 1][GRID_SIZE - 1]);
    }

    // last row on the bottom
    // bottom-left corner element
    NextTemps[PPRows - 1][0] = PPTemps[PPRows - 1][0] + CALC_DTEMP(PPTemps[PPRows - 1][0], left, PPTemps[PPRows - 1][1], PPTemps[PPRows - 2][0], down[0]);

    // middle elements in the last row
    for (int i = 1; i < GRID_SIZE - 1; i++) {
        NextTemps[PPRows - 1][i] = PPTemps[PPRows - 1][i] + CALC_DTEMP(PPTemps[PPRows - 1][i], PPTemps[PPRows - 1][i - 1], PPTemps[PPRows - 1][i + 1], PPTemps[PPRows - 2][i], down[i]);
    }

    // bottom-right corner element
    NextTemps[PPRows - 1][GRID_SIZE - 1] = PPTemps[PPRows - 1][GRID_SIZE - 1] + CALC_DTEMP(PPTemps[PPRows - 1][GRID_SIZE - 1], PPTemps[PPRows - 1][GRID_SIZE - 2], right, PPTemps[PPRows - 2][GRID_SIZE - 1], down[GRID_SIZE - 1]);

    // update the local dataset
    for (int i = 0; i < PPRows; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            PPTemps[i][j] = NextTemps[i][j];
        }
    }
}

void GatherResult(int me) {
    MPI_Status status;

    if (me == BOSS) {
        // Copy the BOSS local temp array data to global temp array
        for (int i = 0; i < PPRows; i++) {
            memcpy(TempData[i], PPTemps[i], GRID_SIZE * sizeof(float));
        }

        // Receive data from other processors
        for (int src = 1; src < NumCpus; src++) {
            int startRow = src * PPRows; // index of the first row of the strip
            int endRow = startRow + PPRows - 1;
            int rowIdx = startRow; // Index of the next row to be received

            while (rowIdx <= endRow) {
                MPI_Recv(TempData[rowIdx], GRID_SIZE, MPI_FLOAT, src, rowIdx, MPI_COMM_WORLD, &status);
                if(DEBUG) fprintf(stderr, "BOSS received %3d from %3d\n", rowIdx, src);
                rowIdx++;
            }
        }

    } else {
        // Send the local temp array PPTemps data back to the BOSS processor row by row
        int startRow = me * PPRows; // Index of the first row of the strip

        for (int i = 0; i < PPRows; i++) {
            int tag = startRow + i; // message tag of MPI_Send
            MPI_Send(PPTemps[i], GRID_SIZE, MPI_FLOAT, BOSS, tag, MPI_COMM_WORLD);
            if(DEBUG) fprintf(stderr, "%3d sent %3d to BOSS\n", me, tag);
        }
    }
}
