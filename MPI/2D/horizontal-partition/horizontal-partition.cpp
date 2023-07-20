#define BOSS 0
#define GRID_SIZE // global 2d array size. Total number of elements = GRID_SIZE X GRID_SIZE

// Heat Diffusion Equation
#define CALC_DTEMP(elem, left, right, up, down) (k_over_rho_c * DT * ((down - 2*elem + up)/(DY * DY) + (left - 2*elem + right)/(DX * DX)))

int NumCpus; // total # of cpus involved
int PPRows;  // per-processor local 2d array size, local number of elements = PPRows X GRID_SIZE
int me;      // which one I am - the rank of a processor

float** PPTemps;   // per-processor local 2d array temperature data
float** NextTemps; // per-processor 2d array to hold next-values
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data
MPI_Status status;

PPRows = GRID_SIZE / NumCpus;

PPTemps = new float* [PPRows];
for (int i = 0; i < PPRows; i++) {
    PPTemps[i] = new float[GRID_SIZE];
}

NextTemps = new float* [PPRows];
for (int i = 0; i < PPRows; i++) {
    NextTemps[i] = new float[GRID_SIZE];
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
            rowIdx++;
        }
    }
}
else {
    // Receive the plate data from the BOSS processor
    int rowIdx = me * PPRows; // Index of the next row to be received

    for (int i = 0; i < PPRows; i++) {
        MPI_Recv(&PPTemps[i][0], GRID_SIZE, MPI_FLOAT, BOSS, rowIdx, MPI_COMM_WORLD, &status);
        rowIdx++;
    }
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
    NextTemps[PPRows - 1][GRID_SIZE - 1] = PPTemps[PPRows - 1][GRID_SIZE - 1] + CALC_DTEMP(PPTemps[PPRows - 1][GRID_SIZE - 1], PPTemps[PPRows - 1][GRID_SIZE - 2], right, PPTemps[PPRows - 2][GRID_SIZE - 1] + down[GRID_SIZE - 1]);

    // update the local dataset
    for (int i = 0; i < PPRows; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            PPTemps[i][j] = NextTemps[i][j];
        }
    }
}
