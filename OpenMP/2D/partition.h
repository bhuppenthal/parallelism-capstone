#define NUM_TIME_STEPS  100                      // number of time steps the simultation runs through

#ifndef SIDE
#define SIDE            8
#endif

const int NUME = SIDE*SIDE;

#ifndef NUMT
#define NUMT            4                       // number of threads to use
#endif

#define NUM_ELEM_PER_THREAD    (NUME/NUMT)      // number of elements in each thread

#ifndef PRINT_ALL_TIME_STEPS
#define PRINT_ALL_TIME_STEPS		true         // set to true to allow all time steps to print
#endif

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

// Partition information and methods
struct partition {
    int row_start;
    int row_end;
    int col_start;
    int col_end;
};

struct partition partitions[NUMT];

void Partition_2D_Array(int Partition_Rows, int Partition_Cols) {
    /*
    Given a partition scheme Partition_Rows and Partition_Cols, determine the rows and columns that
    need to be calculated by each thread.

    A partition scheme of 2 x 6 will result in two rows, six columns.
    ----------------
    |  |  |  |  |  |
    ----------------
    |  |  |  |  |  |
    ----------------
    
    Receives: Partition_Rows: number of divisions across the rows
              Partition_Cols: number of divisions across the columns
    */

    int row_min = SIDE / Partition_Rows;
    int col_min = SIDE / Partition_Cols;

    int row_rem = SIDE % Partition_Rows;
    int col_rem = SIDE % Partition_Cols;

    int row = 0;
    int col = 0;

    for (int i = 0; i < NUMT; i++) {
        partitions[i].row_start = row;
        row += row_min;
        if (i < row_rem)
            row += 1;
        partitions[i].row_end = row;
        row += 1;

        partitions[i].col_start = col;
        col += col_min;
        if (i < col_rem)
            col += 1;
        partitions[i].col_end = row;
        col += 1;
    }
}

void Generate_Partitions() {
    return;
}