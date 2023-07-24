#include <stdio.h>
#include <math.h>
#include <mpi.h>

const float RHO = 8050.;
const float C = 0.466;
const float K = 20.;
float k_over_rho_c = K / (RHO*C);   // units of m^2/sec, cannot be const
// K / (RHO*C) = 5.33x10^-6 m^2/sec

const float DX = 1.0;
const float DT = 1.0;

#define BOSS 0

#define NUMELEMENTS     (8*1024*1024)
#define NUM_TIME_STEPS  4
#define DEBUG           false

float* NextTemps;   // per-processor array to hold computer next-values
int NumCpus;        // total # of cpus involved
int PPSize;         // per-processor local array size
float* PPTemps;     // per-processor local array temperature data
float* TempData;    // the overall NUMELEMENTS-big temperature data

void DoOneTimeStep( int );

int
main( int argc, char *argv[ ] )
{
    MPI_Init( &argc, &argv );

    int me; // which one I am

    MPI_Comm_size( MPI_COMM_WORLD, &NumCpus );
    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    // decide how much data to send to each processor:
    PPSize = NUMELEMENTS / NumCpus;     // assuming it comes out evenly
    PPTemps = new float [PPSize]; // all processors now have this uninitialized Local array
    NextTemps = new float [PPSize]; // all processors now have this uninitialized local array too

    // broadcast the constant:
    MPI_Bcast( (void *)&k_over_rho_c, 1, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

    if( me == BOSS ) // this is the data-creator
    {
        TempData = new float [NUMELEMENTS];
        for( int i = 0; i < NUMELEMENTS; i++ )
            TempData[ i ] = 0.;
        TempData[NUMELEMENTS/2] = 100.;
    }

    MPI_Scatter( TempData, PPSize, MPI_FLOAT, PPTemps, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );

    // all the PPTemps arrays have now been filled
    // do the time steps:

    double time0 = MPI_Wtime( );

    for( int steps = 0; steps < NUM_TIME_STEPS; steps++ )
    {
        // do the computation for one time step:
        DoOneTimeStep( me );

        // ask for all the data:
#ifdef WANT_EACH_TIME_STEPS_DATA
        MPI_Gather( PPTemps, PPSize, MPI_FLOAT, TempData, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );
#endif
    }
#ifndef WANT_EACH_TIME_STEPS_DATA
    MPI_Gather( PPTemps, PPSize, MPI_FLOAT, TempData, PPSize, MPI_FLOAT, BOSS, MPI_COMM_WORLD );
#endif

    double time1 = MPI_Wtime( );

    if( me == BOSS ) {
        double seconds = time1 - time0;
        double performance =
                (double)NUM_TIME_STEPS * (double)NUMELEMENTS / seconds / 1000000.; // mega-elements computed per second
        fprintf( stderr, "%3d, %10d, %8.2lf\n", NumCpus, NUMELEMENTS, performance );
    
    }

    MPI_Finalize( );
    return 0;
}

// read from PerProcessorData[ ], write into NextTemps[ ]
void
DoOneTimeStep( int me ) {
    MPI_Status status;
    // send out the left and right end values:
    // (the tag is from the point of view of the sender)
    if( me != 0 )         // i.e., if i'm not the first group on the left
    {
        // send my PPTemps[0] to me-1 using tag 'L'
        MPI_Send( &PPTemps[0], 1, MPI_FLOAT, me-1, 'L', MPI_COMM_WORLD );
        if( DEBUG ) fprintf( stderr, "%3d sent 'L' to %3d\n", me, me-1 );
    }

    if( me != NumCpus-1 ) // i.e., not the last group on the right
    {
        // send my PPTemps[PPSize-1] to me+1 using tag 'R'
        MPI_Send( &PPTemps[PPSize-1], 1, MPI_FLOAT, me+1, 'R', MPI_COMM_WORLD );
        if( DEBUG ) fprintf( stderr, "%3d sent 'R' to %3d\n", me, me+1 );
    }

    float left = 0.;
    float right = 0.;

    if( me != 0 )     // i.e., if i'm not the first group on the left
    {
        // receive my "left" from me-1 using tag 'R'
        MPI_Recv( &left, 1, MPI_FLOAT, me-1, 'R', MPI_COMM_WORLD, &status );
        if( DEBUG ) fprintf( stderr, "%3d received 'R' from %3d\n", me, me-1 );
    }

    if( me != NumCpus-1 )   // i.e., not the last group on the right
    {
        // receive my "right" from me+1 using tag 'L'
        MPI_Recv( &right, 1, MPI_FLOAT, me+1, 'L', MPI_COMM_WORLD, &status );
        if( DEBUG ) fprintf( stderr, "%3d received 'L' from %3d\n", me, me+1 );
    }

    // first element on the left (0):
    {
        float dtemp = ( k_over_rho_c *
                        ( left - 2.*PPTemps[0] + PPTemps[1] ) / ( DX*DX ) ) * DT;
        NextTemps[0] = PPTemps[0] + dtemp;
    }

    // all the nodes in the middle:
    for( int i = 1; i < PPSize-1; i++ )
    {
        float dtemp = ( k_over_rho_c *
                        ( PPTemps[i-1] - 2.*PPTemps[ i ] + PPTemps[i+1] ) / ( DX*DX ) ) * DT;
        NextTemps[i] = PPTemps[i] + dtemp;
    }

    // last element on the right (PPSize-1):
    {
        float dtemp = ( k_over_rho_c *
                        ( PPTemps[PPSize-2] - 2.*PPTemps[PPSize-1] + right ) / ( DX*DX ) ) * DT;
        NextTemps[PPSize-1] = PPTemps[PPSize-1] + dtemp;
    }

    // update the local dataset:
    for( int i = 0; i < PPSize; i++ )
    {
        PPTemps[ i ] = NextTemps[ i ];
    }
}
