#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mpiio_dmat.h"

void
printMat (double *mat, int size, int rank, int nprocs, char *modif) {
    int i,j;
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    for (j = 0; j < nbR; j++) {
        for (i = 0; i < size; i++) {
            printf ("%s%d %d %lf\n", modif, rank * nbR + j, i, mat[i + j * size]);
        }
    }
}

void
lu_col (int size, double *m, int nprocs, int rank) {
    int i, j, k, k_mod, k_loc, i_incr;
    double tmp, *ik;
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    for (k = 0; k < size - 1; k++ ){
        k_loc = k / nprocs;
        k_mod = k % nprocs;
        ik = (double *) malloc ((size - k - 1) * sizeof (double));

        if(rank == k_mod) {
            for (i = k + 1; i < size; i++){
                m[k_loc * size + i] = m[k_loc * size + i] / m[k_loc * size + k];
                ik[i - k - 1] = m[k_loc * size + i];
            }
        }

        if(rank <= k_mod) {
            k_loc++;
        }

        MPI_Bcast(ik, size - k - 1, MPI_DOUBLE, k_mod, MPI_COMM_WORLD);

        for (j = k_loc; j < nbR; j++){
            tmp = m[j * size + k];
            for (i = k + 1; i < size; i++){
                m[j * size + i] = m[j * size + i] - ik[i - k - 1] * tmp;
            }
        }
        free(ik);
    }
}

void
lu_row (int size, double *m, int nprocs, int rank) {
    int i, j, k, k_mod, k_loc, i_incr;
    double tmp, *akj;
    int nbR = size / nprocs;
    int mod = size % nprocs;
    akj = (double *) malloc (size * sizeof (double));

    if(rank < mod)
        nbR++;

    for (k = 0; k < size - 1; k++ ){
        k_loc = k / nprocs;
        k_mod = k % nprocs;
        if(rank == k_mod) {
            for (i = k; i < size; i++){
                akj[i] = m[k_loc * size + i];
            }
        }
        MPI_Bcast(akj + k, size - k, MPI_DOUBLE, k_mod, MPI_COMM_WORLD);

        if(rank <= k_mod) {
            k_loc++;
        }

        for (i = k_loc; i < nbR; i++) {
            m[i * size + k] /= akj[k];
        }

        for (i = k_loc; i < nbR; i++){
            tmp = m[i * size + k];
            for (j = k + 1; j < size; j++){
                m[i * size + j] = m[i * size + j] - tmp * akj[j];
            }
        }
    }
    free(akj);
}

void
testLoad (int size, int world_size, int world_rank) {
    double *m;
    struct timeval ts, te, t1, t2;

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&ts, 0);
    }

    m = mat_malloc(size, world_rank, world_size);
    mat_read_cyclic(size, world_rank, world_size, m, "a.bin");

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t1, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    lu_row (size, m, world_size, world_rank);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t2, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    mat_write_cyclic(size, world_rank, world_size, m, "lu.bin");
    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&te, 0);
        printf ("%f\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
        printf ("%f\n", (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
        printf ("%f\n", (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
    }

    free (m);
}

int
main (int argc, char **argv) {
    // Initialize the MPI environment
    MPI_Init (NULL, NULL);

    int size;
    if (argc == 2) {
        size = atoi (argv[1]);
    } else {
        size = 16;
    }

    // Get the number of processes
    int world_size;
    MPI_Comm_size (MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);

    testLoad (size, world_size, world_rank);

    // Finalize the MPI environment.
    MPI_Finalize ();
    return 0;
}
