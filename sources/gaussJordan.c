#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define SAVE 1

double *
initMat (int size, int l, int rank) {
    //block distribution of the rows
    double *mat;
    mat = (double *) malloc (l * size * sizeof (double));
    srandom ((unsigned) 233 * rank);
    int i,j;
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
            printf ("%s%d %d %lf\n", modif, rank * nbC + j, i, mat[i + j * size]);
        }
    }
}

double *
loadMat (int size, int l, int rank, char *path) {
    double *mat;
    mat = (double *) malloc (l * size * sizeof (double));
    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at_all (fh, rank * l * size * sizeof (double), mat, l * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    return mat;
}

void
saveMat (double *mat, int size, int l, int rank, char *path) {
    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all (fh, rank * l * size * sizeof (double), mat, l * size, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
}

double *
initVect (int nbC, int rank) {
    double *v;
    int i;
    v = (double *) malloc (nbC * sizeof (double));
    srandom ((unsigned) 21 * rank + 13456);
    for (i = 0; i < nbC; i++) {
        v[i] = 100.0 * rand () / RAND_MAX;
        // printf("v %d %lf\n", i*nbC + rank, v[i]);
    }
    return v;
}

void
printVect (double *v, int nbC, int rank, char *modif) {
    int i;
    for (i = 0; i < nbC; i++) {
        printf ("%s%d %lf\n", modif, nbC * rank + i, v[i]);
    }
}

double *
loadVect (int nbC, int rank, char *path) {
    MPI_File fh;
    MPI_Status status;
    double *v;
    v = (double *) malloc (nbC * sizeof (double));
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at_all (fh, rank * nbC * sizeof (double), v, nbC, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
    return v;
}

void
saveVect (double *v, int nbC, int rank, char *path) {
    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all (fh, rank * nbC * sizeof (double), v, nbC, MPI_DOUBLE, &status);
    MPI_File_close (&fh);
}

void
gaussJordan (int size, double *m, double *v, int nbC, int world_size, int world_rank) {
    int i,j,k;
    double mkk, vk, *akj;
    //TODO : optimize the world_size last steps were the broadcast is not efficient since operations are not performed on all the processes
    for (k = 0; k < size; k++) {
        //kr : local row
        //kw : processus where the data are stored
        int kw = k / nbC;
        int kr = k % nbC;
        if (world_rank == kw) {
            mkk = m[kr * size + k];
            //printf ("mkk %lf kr%d kw%d\n", mkk, kr, kw);
        }
        MPI_Bcast (&mkk, 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

        if (world_rank == kw) {
            //printf ("mkk %lf vk %lf kr%d kw%d\n", mkk, v[kw], kr, kw);
            v[kr] = v[kr] / mkk;
            vk = v[kr];
        }
        MPI_Bcast (&vk, 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

        if (world_rank == kw) {
            for (i = k + 1; i < size; i++) {
                m[kr * size + i] /= mkk;
            }
        }

        akj = (double *) malloc ((size - k - 1) * sizeof (double));

        if (world_rank == kw) {
            for (i = k + 1; i < size; i++) {
                akj[i - k - 1] = m[kr * size + i];
            }
        }
        MPI_Bcast (akj, size - k - 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

        if (world_rank == kw) {
            for (i = 0; i < kr; i++) {
                for (j = k + 1; j < size; j++) {
                    //step 5
                    m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
                }
                v[i] -= m[i * size + k] * vk;
            }
            for (i = kr + 1; i < nbC; i++) {
                for (j = k + 1; j < size; j++) {
                    //step 5
                    m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
                }
                v[i] -= m[i * size + k] * vk;
            }
        } else {
            for (i = 0; i < nbC; i++) {
                for (j = k + 1; j < size; j++) {
                    //step 5
                    m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
                }
                v[i] -= m[i * size + k] * vk;
            }
        }
        free (akj);
    }
}

void
testLoad (int size, int nbC, int world_size, int world_rank) {
    double *v, *m;
    struct timeval ts, te, t1, t2;

    v = initVect (nbC, world_rank);
    //printVect (v, nbC, world_rank, "v - ");
    saveVect (v, nbC, world_rank, "v1.bin");
    free (v);

    m = initMat (size, nbC, world_rank);
    //printMat (m, size, nbC, world_rank, "a - ");
    MPI_Barrier (MPI_COMM_WORLD);
    saveMat (m, size, nbC, world_rank, "m1,1.bin");
    free (m);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&ts, 0);
    }

    m = loadMat (size, nbC, world_rank, "m1,1.bin");
    v = loadVect (nbC, world_rank, "v1.bin");

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t1, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    gaussJordan (size, m, v, nbC, world_size, world_rank);

    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&t2, 0);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    saveVect (v, nbC, world_rank, "v2.bin");
    MPI_Barrier (MPI_COMM_WORLD);
    if (world_rank == 0) {
        gettimeofday (&te, 0);
        printf ("%f\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
        printf ("%f\n", (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
        printf ("%f\n", (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
    }
    //printVect (v, nbC, world_rank, "r - ");
    saveMat (m, size, nbC, world_rank, "m2,2.bin");

    free (m);
    free (v);
}

void
testGen (int size, int nbC, int world_size, int world_rank) {
    double *v, *m;

    v = initVect (nbC, world_rank);
    printVect (v, nbC, world_rank, "v - ");
    saveVect (v, nbC, world_rank, "v1.bin");

    m = initMat (size, nbC, world_rank);
    //printMat (m, size, nbC, world_rank, "a - ");
    MPI_Barrier (MPI_COMM_WORLD);
    saveMat (m, size, nbC, world_rank, "m1,1.bin");

    MPI_Barrier (MPI_COMM_WORLD);
    gaussJordan (size, m, v, nbC, world_size, world_rank);

    MPI_Barrier (MPI_COMM_WORLD);
    saveMat (m, size, nbC, world_rank, "m2,2.bin");
    //printVect (v, nbC, world_rank, "r - ");
    saveVect (v, nbC, world_rank, "v2.bin");

    free (m);
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
