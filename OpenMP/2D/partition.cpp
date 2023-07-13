#include <stdlib.h>

#include "partition.h"

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

    // allocate variable length arrays to hold rows, cols
    int* rows;
    int* cols;
    rows = (int*) malloc(sizeof(int) * Partition_Rows);
    cols = (int*) malloc(sizeof(int) * Partition_Cols);

    for (int i = 0; i < Partition_Rows; i++) {
        rows[i] = row_min;
        if (i < row_rem) rows[i]++;
    }

    for (int j = 0; j < Partition_Cols; j++) {
        cols[j] = col_min;
        if (j < col_rem) cols[j]++;
    }

    int row = 0;
    int thread = 0;
    for (int i = 0; i < Partition_Rows; i++) {
        int col = 0;

        for (int j = 0; j < Partition_Cols; j++) {
            partitions[thread].row_start = row;
            partitions[thread].row_end = row + rows[i] - 1;
            partitions[thread].col_start = col;
            partitions[thread].col_end = col + cols[j] - 1;
            col += cols[j];
            thread++;
        }

        row += rows[i];
    }

    free(rows);
    free(cols);
}