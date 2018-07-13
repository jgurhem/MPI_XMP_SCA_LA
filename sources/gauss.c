#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>

#define SAVE 1
// SAVE==0 => no save

double *
initMat (int size, int l, int rank) {
    double *mat;
    int i,j;
    mat = (double *) malloc (l * size * sizeof (double));
    srandom ((unsigned) 233 * rank);
    for (j = 0; j < l; j++) {
        for (i = 0; i < size; i++) {
            mat[i + j * size] = 100.0 * rand () / RAND_MAX;;
            // printf("a %d %d %lf\n", i, rank+l*j, mat[i*size+j]); 
        }
    }
    return mat;
}

void
printMat (double *mat, int size, int nbC, int rank, char *modif) {
    int i,j;
    for (j = 0; j < nbC; j++) {
        for (i = 0; i < size; i++) {
            // cyclic distribution of the columns of the matrix
            printf ("%s%d %d %lf\n", modif, i, rank + j * nbC, mat[i + j * size]);
        }
    }
}

double *
loadMat (int size, int l, int rank, char *path) {
    int world_size;
    double *mat;
    mat = (double *) malloc (l * size * sizeof (double));
    MPI_Comm_size (MPI_COMM_WORLD, &world_size);

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[2] = { size, size };
    int array_distrib[2] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_NONE };
    int array_dargs[2] = { 1, size };
    int array_psizes[2] = { world_size, 1 };
    MPI_Type_create_darray (world_size /* size */ , rank, 2 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_read_all (fh, mat, l * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
    return mat;
}

void
saveMat (double *mat, int size, int l, int rank, char *path) {
    int world_size;
    MPI_Comm_size (MPI_COMM_WORLD, &world_size);

    MPI_File fh;
    MPI_Status status;
    MPI_Datatype darray;
    int array_size[2] = { size, size };
    int array_distrib[2] = { MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_NONE };
    int array_dargs[2] = { 1, size };
    int array_psizes[2] = { world_size, 1 };
    MPI_Type_create_darray (world_size /* size */ , rank, 2 /* dims */ , array_size, array_distrib, array_dargs, array_psizes, MPI_ORDER_C, MPI_DOUBLE, &darray);
    MPI_Type_commit (&darray);

    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view (fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_write_all (fh, mat, l * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    MPI_Type_free (&darray);
}

double *
initVect (int size) {
    double *v;
    int i;
    v = (double *) malloc (size * sizeof (double));
    srandom ((unsigned) 23);
    for (i = 0; i < size; i++) {
        v[i] = 100.0 * rand () / RAND_MAX;
        // printf("v %d %lf\n", i, v[i]);
    }
    return v;
}

void
printVect (double *v, int size, char *modif) {
    int i;
    for (i = 0; i < size; i++) {
        printf ("%s%d %lf\n", modif, i, v[i]);
    }
}

double *
loadVect (int size, char *path) {
    int f;
    double *v;
    v = (double *) malloc (size * sizeof (double));
    f = open (path, O_RDONLY);
    read (f, v, size * sizeof (double));
    close (f);
    return v;
}

void
saveVect (double *v, int size, char *path) {
    int f;
    f = open (path, O_WRONLY | O_CREAT | O_TRUNC, 0664);
    write (f, v, size * sizeof (double));
    close (f);
}

void
gaussEliminationMPI (int size, double *m, double *v, int nbC, int world_size, int world_rank) {

    MPI_Status Stat;
    int i,j,k;
    double mkk, *aik;
    //TODO : optimize the world_size last steps were the broadcast is not efficient since operations are not performed on all the processes
    for (k = 0; k < size; k++) {
        //kw : local column
        //kr : processus where the data are stored for the cyclic distribution
        int kw = k / world_size;
        int kr = k % world_size;

        if (world_rank == kr) {
            mkk = m[kw * size + k];
        }
        MPI_Bcast (&mkk, 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);

        if (world_rank == 0) {
            v[k] = v[k] / mkk;
        }

        if (world_rank > kr) {
            for (i = kw; i < nbC; i++) {
                m[i * size + k] /= mkk;
            }
        } else {
            for (i = kw + 1; i < nbC; i++) {
                m[i * size + k] /= mkk;
            }
        }

        aik = (double *) malloc ((size - k - 1) * sizeof (double));

        if (world_rank == kr) {
            for (i = k + 1; i < size; i++) {
                aik[i - k - 1] = m[kw * size + i];
            }
        }
        MPI_Bcast (aik, size - k - 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);

        if (world_rank > kr) {
            for (j = kw; j < nbC; j++) {
                for (i = k + 1; i < size; i++) {
                    //step 5
                    m[j * size + i] -= aik[i - k - 1] * m[j * size + k];
                }
            }
        } else {
            for (j = kw + 1; j < nbC; j++) {
                for (i = k + 1; i < size; i++) {
                    //step 5
                    m[j * size + i] -= aik[i - k - 1] * m[j * size + k];
                }
            }
        }
        if (world_rank == 0) {
            //step 4
            for (i = k + 1; i < size; i++) {
                v[i] -= aik[i - k - 1] * v[k];
            }
        }

        free (aik);
    }

    for (k = size - 1; k > 0; k--) {
        //k = 7;
        //kw : local column
        //kr : processus where the data are stored for the cyclic distribution
        int kw = k / world_size;
        int kr = k % world_size;

        if (kr == 0 && world_rank == 0) {
            for (i = 0; i < k; i++) {
                v[i] -= m[kw * size + i] * v[k];
            }
            continue;
        }

        if (world_rank == kr) {
            aik = (double *) malloc (k * sizeof (double));
            for (i = 0; i < k; i++) {
                aik[i] = m[kw * size + i];
            }
            MPI_Send (aik, k, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
            free (aik);
        }

        if (world_rank == 0) {
            aik = (double *) malloc (k * sizeof (double));
            MPI_Recv (aik, k, MPI_DOUBLE, kr, k, MPI_COMM_WORLD, &Stat);
            for (i = 0; i < k; i++) {
                v[i] -= aik[i] * v[k];
            }
            free (aik);
        }
    }

}

void
testGen (int size, int nbC, int world_size, int world_rank) {
    double *v, *m;

    if (world_rank == 0) {
        v = initVect (size);
        //printVect (v, size,"v - ");
        if (SAVE) {
            saveVect (v, size, "v1.bin");
        }
    }

    m = initMat (size, nbC, world_rank);
    //printMat (m, size,nbC, world_rank, "a - ");
    if (SAVE) {
        saveMat (m, size, nbC, world_rank, "m1,1.bin");
    }

    gaussEliminationMPI (size, m, v, nbC, world_size, world_rank);

    //printMat (m, size,nbC, world_rank, "m - ");
    if (SAVE) {
        saveMat (m, size, nbC, world_rank, "m2,2.bin");

        if (world_rank == 0) {
            //printVect (v,size, "r - ");
            saveVect (v, size, "v2.bin");
        }
    }

    free (m);
    if (world_rank == 0)
        free (v);
}

void
testLoad (int size, int nbC, int world_size, int world_rank) {
    double *v, *m;
    struct timeval ts, te, t1, t2;

    if (world_rank == 0) {
        v = initVect (size);
        //printVect (v, size,"v - ");
        saveVect (v, size, "v1.bin");
        free (v);
    }

    m = initMat (size, nbC, world_rank);
    saveMat (m, size, nbC, world_rank, "m1,1.bin");
    free (m);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&ts, 0);
    }

    m = loadMat (size, nbC, world_rank, "m1,1.bin");
    if (world_rank == 0) {
        v = loadVect (size, "v1.bin");
    }

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t1, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);

    gaussEliminationMPI (size, m, v, nbC, world_size, world_rank);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t2, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);

    //printMat (m, size,nbC, world_rank, "m - ");
    saveMat (m, size, nbC, world_rank, "m2,2.bin");

    if (world_rank == 0) {
        //printVect (v,size, "r - ");
        saveVect (v, size, "v2.bin");
        gettimeofday (&te, 0);
        printf ("%f\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
        printf ("%f\n", (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
        printf ("%f\n", (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
    }

    free (m);
    if (world_rank == 0)
        free (v);
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

    // number of columns of the matrix distributed per processes
    int nbC = size / world_size;
    assert (nbC * world_size == size);

    testLoad (size, nbC, world_size, world_rank);

    // Finalize the MPI environment.
    MPI_Finalize ();
    return 0;
}
