#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <iostream>

const float RHO = 8050.;
const float C = 0.466;
const float K = 20.;
float k_over_rho_c = K / (RHO*C);   // units of m^2/sec, cannot be const
// K / (RHO*C) = 5.33x10^-6 m^2/sec

const float DX = 1.0;
const float DY = 1.0;
const float DT = 1.0;

#define BOSS 0

#define NUM_TIME_STEPS  100
#define DEBUG           false
// #define NUMELEMENTS     (8*1024*1024)
// #define NUMELEMENTS     8 * 8
// #define SIDE            sqrt(NUMELEMENTS)
#define SIDE            8
#define NUMELEMENTS     SIDE * SIDE

// Heat Diffusion Equation
#define CALC_DTEMP(elem, left, right, up, down) (k_over_rho_c * DT * ((down - 2*elem + up)/(DY * DY) + (left - 2*elem + right)/(DX * DX)));

#define WANT_EACH_TIME_STEPS_DATA       true

float* NextTemps;   // per-processor array to hold computer next-values
int NumCpus;        // total # of cpus involved
int PPSize;         // per-processor local array size
int PPSide;         // side of each processors sub-square in the grid of squares partition
float* PPTemps;     // per-processor local array temperature data
float** TempData;    // the overall NUMELEMENTS-big temperature data
float* OneDimPartitionTemp; // 1D array to copy into

void DoOneTimeStep( int );

int
main( int argc, char *argv[ ] )
{
    // fprintf(stderr, "!!! started program\n");

    MPI_Init( &argc, &argv );

    int me; // which one I am

    MPI_Comm_size( MPI_COMM_WORLD, &NumCpus );
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    // decide how much data to send to each processor:
    PPSize = NUMELEMENTS / NumCpus;     // assuming it comes out evenly
    PPSide = sqrt(PPSize);              // also assuming it comes out evenly
    // PPTemps = new float [sqrt(PPSize)] [sqrt(PPSize)]; // all processors now have this uninitialized Local array
    // NextTemps = new float [sqrt(PPSize)] [sqrt(PPSize)]; // all processors now have this uninitialized local array too
    PPTemps = new float [PPSize]; // all processors now have this uninitialized Local array
    NextTemps = new float [PPSize]; // all processors now have this uninitialized local array too

    // fprintf(stderr, "!!! passed NextTemps allocation\n");

    // broadcast the constant:
    MPI_Bcast( (void *)&k_over_rho_c, 1, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

    if( me == BOSS ) // this is the data-creator
    {
        // create 2D array with 100 temp in the middle element and 0 everywhere else
        // TempData = new float [SIDE][SIDE];
        TempData = new float* [SIDE];

        // fprintf(stderr, "%016x\n", TempData);
        // fprintf(stderr, "%016x\n", &TempData);
        // fprintf(stderr, "%d\n", TempData[0]);
        // fprintf(stderr, "%016x\n", &(&TempData));
        // fprintf(stderr, "!!! passed TempData allocation\n");

        for( int i = 0; i < SIDE; i++ )
        {
            TempData[i] = new float [SIDE];
            for( int j = 0; j < SIDE; j++ )
            {
                TempData[ i ][ j ] = 0.;
                // TempData[ i ][ j ] = i * SIDE + j;
                // fprintf(stderr, "!!! passed inner for loop filling data\n");
            }
        }
        // TempData[SIDE/2][SIDE/2] = 100.;
        TempData[SIDE/2 - 1][SIDE/2 - 1] = 100.;

        for (int i = 0; i < SIDE; i++) {
            for (int j = 0; j < SIDE; j++) {
                fprintf(stderr, " %.2f ", TempData[i][j]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n\n");

        // fprintf(stderr, "!!! passed TempData filling data\n");

        // copy each sub-square of the grid of squares partition into a 1D array,
        // in which the sections of the array sent to each processor are sequential and contiguous
        OneDimPartitionTemp = new float [NUMELEMENTS];

        int k = 0;
        for (int starting_row = 0; starting_row < SIDE; starting_row += PPSide)
        {
            for (int starting_col = 0; starting_col < SIDE; starting_col += PPSide)
            {
                for (int i = starting_row; i < starting_row + PPSide; i++)
                {
                    for (int j = starting_col; j < starting_col + PPSide; j++)
                    {
                        OneDimPartitionTemp[k] = TempData[i][j];
                        k++;
                    }
                }
            }
        }

        // int starting_row = 0;
        // int starting_column = 0;
        // for (int i = 0; i < NumCpus; i++)
        // {
        //     for (int r = 0; r < PPSide; r++)
        //     {
        //         for (int c = 0; c < PPSide; c++)
        //         {
        //             int one_d_array_index = i * PPSize + r * PPSide + c;    
        //             OneDimPartitionTemp[one_d_array_index] = TempData[starting_row + r][starting_column + c];
        //         }
        //     }

        //     if ((i + 1) % (SIDE / PPSide) == 0)
        //     {
        //         starting_row = starting_row + PPSide;
        //     }
        // }

        // for (int i=0; i<NUMELEMENTS; i++)
        // {
        //     fprintf(stderr, "%0.2f ", OneDimPartitionTemp[i]);
        // }
        // fprintf(stderr, "\n");
    }

    if (me == BOSS) {
        // fprintf(stderr, "!!! passed TempData copy to 1D\n");
    }
    
    // MPI_Scatter( TempData, PPSize, MPI_FLOAT, PPTemps, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );
    MPI_Scatter( OneDimPartitionTemp, PPSize, MPI_FLOAT, PPTemps, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

    if (me == BOSS) {
        // fprintf(stderr, "!!! passed Scatter\n");
    }

    // all the PPTemps arrays have now been filled
    // do the time steps:

    double time0 = MPI_Wtime( );

    for( int steps = 0; steps < NUM_TIME_STEPS; steps++ )
    {
        // do the computation for one time step:
        DoOneTimeStep( me );

        // ask for all the data:
#ifdef WANT_EACH_TIME_STEPS_DATA
        MPI_Gather( PPTemps, PPSize, MPI_FLOAT, OneDimPartitionTemp, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

        if (me == BOSS) {
            // fprintf(stderr, "!!! passed Gather Each Time Step\n");
        }

        if (me == BOSS) {
            // convert OneDimPartitionTemp back into 2D array
            int k = 0;
            for (int starting_row = 0; starting_row < SIDE; starting_row += PPSide)
            {
                for (int starting_col = 0; starting_col < SIDE; starting_col += PPSide)
                {
                    for (int i = starting_row; i < starting_row + PPSide; i++)
                    {
                        for (int j = starting_col; j < starting_col + PPSide; j++)
                        {
                            TempData[i][j] = OneDimPartitionTemp[k];
                            k++;
                        }
                    }
                }
            }

            // int k = 0;      // where to get the temp from 1D array
            // for( int i = 0; i < SIDE; i++ )
            // {
            //     for( int j = 0; j < SIDE; j++ )
            //     {
            //         {
            //             TempData[i][j] = OneDimPartitionTemp[k];
            //             k++;
            //         }
            //     }
            // }
        }

        if (me == BOSS) {
            // fprintf(stderr, "!!! passed covert 1D back to 2D\n");
        }

        // only the Boss prints data
        if( me == BOSS ) {
            // // print temperature data for verification
            // fprintf(stderr, "Time step: %i ", steps);

            // for (int i = 0; i < NUMELEMENTS; i++) {
            //     fprintf(stderr, " %.2f ", TempData[i]);
            // }
            // fprintf(stderr, "\n");

            // // Print out TempData to stderr
            // fprintf(stdout, "Initial Matrix\n");
            // for (int i = 0; i < SIDE; i++) {
            //     for (int j = 0; j < SIDE; j++) {
            //         std::cout << std::fixed << std::setprecision(2) << TempData[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            for (int i = 0; i < SIDE; i++) {
                for (int j = 0; j < SIDE; j++) {
                    fprintf(stderr, " %.2f ", TempData[i][j]);
                }
                fprintf(stderr, "\n");
            }
            fprintf(stderr, "\n\n");
        }

        if (me == BOSS) {
            // fprintf(stderr, "!!! passed print verification temps\n");
        }


#endif
    }
#ifndef WANT_EACH_TIME_STEPS_DATA
    MPI_Gather( PPTemps, PPSize, MPI_FLOAT, OneDimPartitionTemp, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

    if (me == BOSS) {
        // convert OneDimPartitionTemp back into 2D array
        int k = 0;
        for (int starting_row = 0; starting_row < SIDE; starting_row += PPSide)
        {
            for (int starting_col = 0; starting_col < SIDE; starting_col += PPSide)
            {
                for (int i = starting_row; i < starting_row + PPSide; i++)
                {
                    for (int j = starting_col; j < starting_col + PPSide; j++)
                    {
                        TempData[i][j] = OneDimPartitionTemp[k];
                        k++;
                    }
                }
            }
        }
    }

    if (me == BOSS) {
        // fprintf(stderr, "!!! passed no time steps copy back into 2D\n");
    }

#endif

    if (me == BOSS) {
        // fprintf(stderr, "!!! passed endif\n");
    }

    double time1 = MPI_Wtime( );

    if( me == BOSS ) {
        double seconds = time1 - time0;
        double performance =
                (double)NUM_TIME_STEPS * (double)NUMELEMENTS / seconds / 1000000.; // mega-elements computed per second
        fprintf( stderr, "%3d, %10d, %8.2lf\n", NumCpus, NUMELEMENTS, performance );
    
    }

    if (me == BOSS) {
        // fprintf(stderr, "!!! passed performance calc\n");
    }

    MPI_Finalize( );
    return 0;
}

// read from PerProcessorData[ ], write into NextTemps[ ]
void
DoOneTimeStep( int me ) {

    // for (int i = 0; i < PPSize; i++)
    // {
    //     fprintf(stderr, " %.2f ", PPTemps[i]);
    // }
    // fprintf(stderr, "\n");

    // fprintf(stderr, "!!! rank %d entered DoOneTimeStep\n", me);

    MPI_Status status;
    // send out the left and right end values:
    // send out the top and bottom end values:
    // rank of processors is starting from top left corner and counting left to right across the rows
    // (the tag is from the point of view of the sender)
    // if( me != 0 )         // i.e., if i'm not the first group on the left
    if( me % (SIDE / PPSide) != 0)  // i.e., if i'm not a partition on the leftmost column of the object
    {
        // send my leftmost column of elements to me-1 using tag 'L'
        // MPI_Send( &PPTemps[0], 1, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD );
        // if( DEBUG ) fprintf( stderr, "%3d sent 'L' to %3d\n", me, me-1 );
        for( int i = 0; i < PPSize; i++ )   // loop over elements of partition
        {
            if( i % PPSide == 0)    // if element is on the leftmost column of the partition
            {

                // fprintf(stderr, "!!! PPTemps: %016x\n", PPTemps);
                // fprintf(stderr, "!!! &PPTemps[i]: %016x\n", &PPTemps[i]);

                MPI_Send( &PPTemps[i], 1, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD );
            }
        }
    }

    // if( me != NumCpus-1 ) // i.e., not the last group on the right
    if( me % (SIDE / PPSide) != (SIDE / PPSide) - 1 ) // i.e., if i'm not a partition on the rightmost column of the object
    {
        // send my rightmost column of elements to me+1 using tag 'R'
        // MPI_Send( &PPTemps[PPSize-1], 1, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD );
        // if( DEBUG ) fprintf( stderr, "%3d sent 'R' to %3d\n", me, me+1 );
        for( int i = 0; i < PPSize; i++ )   // loop over elements of partition
        {
            if( i % PPSide == PPSide - 1 )   // if element is on the rightmost column of the partition
            {
                MPI_Send( &PPTemps[i], 1, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD );
            }
        }
    }

    if( me >= (SIDE / PPSide) )     // i.e., if i'm not a partition on the top row of the object
    {
        // send my top column of elements to me-(SIDE / PPSide) using tag 'T'
        // MPI_Send( &PPTemps[0], 1, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD );
        // if( DEBUG ) fprintf( stderr, "%3d sent 'L' to %3d\n", me, me-1 );
        for( int i = 0; i < PPSize; i++ )   // loop over elements of partition
        {
            if( i < PPSide )   // if element is on the top row of the partition
            {
                MPI_Send( &PPTemps[i], 1, MPI_FLOAT, me-(SIDE / PPSide), 'T', MPI_COMM_WORLD ); // me-(SIDE / PPSide) is partition directly above current partition
            }
        }
    }

    if( me < (NumCpus - (SIDE / PPSide))) // i.e., if i'm not a partition on the bottom row of the object
    {
        // send my bottom column of elements to me+(SIDE / PPSide) using tag 'B'
        // MPI_Send( &PPTemps[PPSize-1], 1, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD );
        // if( DEBUG ) fprintf( stderr, "%3d sent 'R' to %3d\n", me, me+1 );
        for( int i = 0; i < PPSize; i++ )   // loop over elements of partition
        {
            if( i >= PPSize - PPSide )   // if element is on the bottom row of the partition
            {
                MPI_Send( &PPTemps[i], 1, MPI_FLOAT, me+(SIDE / PPSide), 'B', MPI_COMM_WORLD ); // me+(SIDE / PPSide) is partition directly below current partition
            }
        }
    }

    // fprintf(stderr, "!!! Rank %d passed send calls\n", me);

    // float left = 0.;
    // float right = 0.;

    float* left = new float [PPSide];
    float* right = new float [PPSide];
    float* up = new float [PPSide];
    float* down = new float [PPSide];

    for( int i = 0; i < PPSide; i++ ) 
    {
        left[i] = 0.;
        right[i] = 0.;
        up[i] = 0.;
        down[i] = 0.;
    }

    // if( me != 0 )     // i.e., if i'm not the first group on the left
    if( me % (SIDE / PPSide) != 0)  // i.e., if i'm not a partition on the leftmost column of the object
    {
        // receive my "left" from me-1 using tag 'R'
        // MPI_Recv( &left, 1, MPI_FLOAT, me-1, 'R', MPI_COMM_WORLD, &status );
        // if( DEBUG ) fprintf( stderr, "%3d received 'R' from %3d\n", me, me-1 );
        for( int i = 0; i < PPSide; i++ )
        {

            // fprintf(stderr, "!!! left: %016x\n", left);
            // fprintf(stderr, "!!! &left[i]: %016x\n", &left[i]);

            MPI_Recv( &left[i], 1, MPI_FLOAT, me-1, 'R', MPI_COMM_WORLD, &status );
        }
    }

    // if( me != NumCpus-1 )   // i.e., not the last group on the right
    if( me % (SIDE / PPSide) != (SIDE / PPSide) - 1 ) // i.e., if i'm not a partition on the rightmost column of the object
    {
        // receive my "right" from me+1 using tag 'L'
        // MPI_Recv( &right, 1, MPI_FLOAT, me+1, 'L', MPI_COMM_WORLD, &status );
        // if( DEBUG ) fprintf( stderr, "%3d received 'L' from %3d\n", me, me+1 );
        for( int i = 0; i < PPSide; i++ )
        {
            MPI_Recv( &right[i], 1, MPI_FLOAT, me+1, 'L', MPI_COMM_WORLD, &status );
        }
    }

    if( me >= (SIDE / PPSide) ) // i.e., if i'm not a partition on the top column of the object
    {
        // receive my "up" from me-(SIDE / PPSide) using tag 'B'
        // me-(SIDE / PPSide) is partition directly above current partition
        // MPI_Recv( &right, 1, MPI_FLOAT, me+1, 'L', MPI_COMM_WORLD, &status );
        // if( DEBUG ) fprintf( stderr, "%3d received 'L' from %3d\n", me, me+1 );
        for( int i = 0; i < PPSide; i++ )
        {
            MPI_Recv( &up[i], 1, MPI_FLOAT, me-(SIDE / PPSide), 'B', MPI_COMM_WORLD, &status );
        }
    }

    if( me < (NumCpus - (SIDE / PPSide))) // i.e., if i'm not a partition on the bottom column of the object
    {
        // receive my "right" from me+1 using tag 'L'
        // me+(SIDE / PPSide) is partition directly below current partition
        // MPI_Recv( &right, 1, MPI_FLOAT, me+1, 'L', MPI_COMM_WORLD, &status );
        // if( DEBUG ) fprintf( stderr, "%3d received 'L' from %3d\n", me, me+1 );
        for( int i = 0; i < PPSide; i++ )
        {
            MPI_Recv( &down[i], 1, MPI_FLOAT, me+(SIDE / PPSide), 'T', MPI_COMM_WORLD, &status );
        }
    
        // if (me == 1)
        // {
        //     fprintf(stderr, "rank 1 down[]\n");
        //     for (int i=0; i< PPSize; i++)
        //     {
        //         fprintf(stderr, "%0.2f ", down[i]);
        //     }
        //     fprintf(stderr, "\n");
        // }
    }

    // fprintf(stderr, "!!! rank %d passed receives\n", me);

    float left_elem;
    float right_elem;
    float up_elem;
    float down_elem;

    // loop through rows of the partition (already contiguous in 1D array for this partition)
    // for (int i = first_row; i <= last_row; i++)
    for (int row = 0; row < PPSide; row++) {
        // leftmost element

        // leftmost element must access column to the left of partition
        // left array is already ordered by row num
        left_elem = left[row];
        // element to the right of leftmost element always accesses element within partition (unless partition is only one element)
        // access element to the right of leftmost element
        right_elem = PPTemps[row * PPSide + 1];

        // if top row of partition, access row above partition
        // up array is already ordered by row num
        up_elem = up[0];
        // if not the top row of the partition
        if (row != 0) {
            // access element above leftmost element within partition
            up_elem = PPTemps[(row - 1) * PPSide];
        }

        // if bottom row of partition, access row below partition
        // down array is already ordered by row num
        down_elem = down[0];
        // if not the bottom row of the partition
        if (row != PPSide - 1) {
            // access element below leftmost element within partition
            down_elem = PPTemps[(row + 1) * PPSide];
        }
        
        // if (me ==1 && row == 3){
        //     fprintf(stderr, "!!! %0.2f ", down_elem);
        // }

        // NextTemps[i][0] = PPTemps[i][0] + CALC_DTEMP(PPTemps[i][0], left[i], right[i], up[i], down[i]);
        NextTemps[row * PPSide] = PPTemps[row * PPSide] + CALC_DTEMP(PPTemps[row * PPSide], left_elem, right_elem, up_elem, down_elem);

        // fprintf(stderr, "!!! rank %d passed leftmost elements\n", me);

        // middle elements
        for (int column = 1; column < PPSide - 1; column++){
            // Arr[i][j] = Arr[ (i*col) + j ]
            left_elem = PPTemps[(row * PPSide) + column - 1];
            right_elem = PPTemps[(row * PPSide) + column + 1];

            // if top row of partition, access row above partition
            // up array is already ordered by column num
            up_elem = up[column];
            // if not the top row of the partition
            if (row != 0) {
                // access element above leftmost element within partition
                up_elem = PPTemps[(row - 1) * PPSide + column];
            }

            // if bottom row of partition, access row below partition
            // down array is already ordered by column num
            down_elem = down[column];
            // if not the bottom row of the partition
            if (row != PPSide - 1) {
                // access element below leftmost element within partition
                down_elem = PPTemps[(row + 1) * PPSide + column];
            }

            NextTemps[(row * PPSide) + column] = PPTemps[(row * PPSide) + column] + CALC_DTEMP(PPTemps[(row * PPSide) + column], left_elem, right_elem, up_elem, down_elem);
        }

        // fprintf(stderr, "!!! rank %d passed middle elements\n", me);

        // rightmost element
        
        // element to the left of rightmost element always accesses element within partition (unless partition is only one element)
        // left array is already ordered by row num
        left_elem = PPTemps[row * PPSide - 1];
        // rightmost element must access column to the right of partition
        // access element to the right of rightmost element
        right_elem = right[row];

        // if top row of partition, access row above partition
        // up array is already ordered by row num
        up_elem = up[PPSide - 1];
        // if not the top row of the partition
        if (row != 0) {
            // access element above rightmost element within partition
            up_elem = PPTemps[(row - 1) * PPSide + PPSide - 1];
        }

        // if bottom row of partition, access row below partition
        // down array is already ordered by row num
        down_elem = down[PPSide - 1];
        // if not the bottom row of the partition
        if (row != PPSide - 1) {
            // access element below rightmost element within partition
            down_elem = PPTemps[(row + 1) * PPSide + PPSide - 1];
        }

        NextTemps[row * PPSide + PPSide - 1] = PPTemps[row * PPSide + PPSide - 1] + CALC_DTEMP(PPTemps[row * PPSide + PPSide - 1], left_elem, right_elem, up_elem, down_elem);
    }

    // fprintf(stderr, "!!! rank %d passed rightmost elements\n", me);

    // update the local dataset:
    for( int i = 0; i < PPSize; i++ )
    {
        PPTemps[ i ] = NextTemps[ i ];
    }

    // fprintf(stderr, "!!! rank %d passed update local data\n", me);
}
