#include <mpi.h>
#include <stdlib.h>

double *mat_malloc (int size, int rank, int nprocs) {
    int nbR = size / nprocs;
    int mod = size % nprocs;
    double *mat;
    if(rank < mod)
        nbR++;
    mat = (double *) malloc (nbR * size * sizeof (double));
    return mat;
}

void mat_read_cyclic(int size, int rank, int nprocs, void *mat, char *path) {
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

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_read_all (fh, mat, nbR * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
}


void mat_write_cyclic(int size, int rank, int nprocs, void *mat, char *path) {
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

void vect_read_cyclic(int size, int rank, int nprocs, void *vect, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[1] = { size };
    int array_distrib[1] = { MPI_DISTRIBUTE_CYCLIC };
    int array_dargs[1] = { 1 };
    int array_psizes[1] = { nprocs };
    MPI_Type_create_darray (nprocs /* size */ , rank, 1 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_read_all (fh, vect, nbR, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
}


void vect_write_cyclic(int size, int rank, int nprocs, void *vect, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[1] = { size };
    int array_distrib[1] = { MPI_DISTRIBUTE_CYCLIC };
    int array_dargs[1] = { 1 };
    int array_psizes[1] = { nprocs };
    MPI_Type_create_darray (nprocs /* size */ , rank, 1 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_write_all (fh, vect, nbR, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
}

void mat_read_block (int size, int rank, int nprocs, void *mat, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(rank < mod) {
        MPI_File_read_at_all (fh, rank * nbR * size * sizeof (double), mat, nbR * size, MPI_DOUBLE, &status);
    } else {
        MPI_File_read_at_all (fh, (rank * nbR + mod) * size * sizeof (double), mat, nbR * size, MPI_DOUBLE, &status);
    }
    MPI_File_close (&fh);
}

void mat_write_block (int size, int rank, int nprocs, void *mat, char *path) {
    int nbR = size / nprocs;
    int mod = size % nprocs;

    if(rank < mod)
        nbR++;

    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_WORLD, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if(rank < mod)
        MPI_File_write_at_all (fh, rank * nbR * size * sizeof (double), mat, nbR * size, MPI_DOUBLE, &status);
    else
        MPI_File_write_at_all (fh, (rank * nbR + mod) * size * sizeof (double), mat, nbR * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
}
