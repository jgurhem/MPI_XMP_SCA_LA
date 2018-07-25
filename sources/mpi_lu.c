#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

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

double *
loadMat (int size, int rank, int nprocs, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;
    double *mat;

    if(rank < mod)
        nbR++;

    mat = (double *) malloc (nbR * size * sizeof (double));

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[2] = { size, size };
    int array_distrib[2] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_NONE };
    int array_dargs[2] = { 1, size };
    int array_psizes[2] = { nprocs, 1 };
    MPI_Type_create_darray (nprocs /* size */ , rank, 2 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_read_all (fh, mat, nbR * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
    return mat;
}

void
saveMat (double *mat, int size, int rank, int nprocs, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[2] = { size, size };
    int array_distrib[2] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_NONE };
    int array_dargs[2] = { 1, size };
    int array_psizes[2] = { nprocs, 1 };
    MPI_Type_create_darray (nprocs /* size */ , rank, 2 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_write_all (fh, mat, nbR * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
}

void
lu (int size, double *m, int nprocs, int rank) {
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
testLoad (int size, int world_size, int world_rank) {
    double *m;
    struct timeval ts, te, t1, t2;

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&ts, 0);
    }

    m = loadMat (size, world_rank, world_size, "a.bin");

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t1, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    lu (size, m, world_size, world_rank);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t2, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    saveMat (m, size, world_rank, world_size, "lu.bin");
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
