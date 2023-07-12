#define BOSS 0
#define GRID_SIZE // global 2d array size. Total number of elements = GRID_SIZE X GRID_SIZE

int NumCpus; // total # of cpus involved
int PPRows;  // per-processor local 2d array size, local number of elements = PPRows X GRID_SIZE
int me;      // which one I am - the rank of a processor

float** PPTemps;   // per-processor local 2d array temperature data
float** TempData;  // the overall 2d array (GRID_SIZE X GRID_SIZE)-big temperature data
MPI_Status status;

PPRows = GRID_SIZE / NumCpus;

PPTemps = new float* [PPRows];
for (int i = 0; i < PPRows; i++) {
    PPTemps[i] = new float[GRID_SIZE];
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
