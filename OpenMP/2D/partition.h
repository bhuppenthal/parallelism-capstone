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
const float k_over_rho_c = K/(RHO*C);                 // units of m^2/s NOTE: Cannot be a const (true for OpenMP too?)
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

extern struct partition partitions[NUMT];

// Function prototypes.
void Partition_2D_Array(int Partition_Rows, int Partition_Cols);