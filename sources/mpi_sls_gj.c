#include "mpiio_dmat.h"
#include "parse_args.h"
#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void printMat(double *mat, int size, int nbC, int rank, char *modif) {
  int i, j;
  for (j = 0; j < nbC; j++) {
    for (i = 0; i < size; i++) {
      printf("%s%d %d %lf\n", modif, rank * nbC + j, i, mat[i + j * size]);
    }
  }
}

double *initVect(int nbC, int rank) {
  double *v;
  int i;
  v = (double *)malloc(nbC * sizeof(double));
  srandom((unsigned)21 * rank + 13456);
  for (i = 0; i < nbC; i++) {
    v[i] = 100.0 * rand() / RAND_MAX;
    // printf("v %d %lf\n", i*nbC + rank, v[i]);
  }
  return v;
}

void printVect(double *v, int nbC, int rank, char *modif) {
  int i;
  for (i = 0; i < nbC; i++) {
    printf("%s%d %lf\n", modif, nbC * rank + i, v[i]);
  }
}

double *loadVect(int nbC, int rank, char *path) {
  MPI_File fh;
  MPI_Status status;
  double *v;
  v = (double *)malloc(nbC * sizeof(double));
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_read_at_all(fh, rank * nbC * sizeof(double), v, nbC, MPI_DOUBLE,
                       &status);
  MPI_File_close(&fh);
  return v;
}

void saveVect(double *v, int nbC, int rank, char *path) {
  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  MPI_File_write_at_all(fh, rank * nbC * sizeof(double), v, nbC, MPI_DOUBLE,
                        &status);
  MPI_File_close(&fh);
}

void gaussJordan(int size, double *m, double *v, int nbC, int world_size,
                 int world_rank) {
  int i, j, k;
  double mkk, vk, *akj;
  // TODO : optimize the world_size last steps were the broadcast is not
  // efficient since operations are not performed on all the processes
  for (k = 0; k < size; k++) {
    // kr : local row
    // kw : processus where the data are stored
    int kw = k / nbC;
    int kr = k % nbC;
    if (world_rank == kw) {
      mkk = m[kr * size + k];
      // printf ("mkk %lf kr%d kw%d\n", mkk, kr, kw);
    }
    MPI_Bcast(&mkk, 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

    if (world_rank == kw) {
      // printf ("mkk %lf vk %lf kr%d kw%d\n", mkk, v[kw], kr, kw);
      v[kr] = v[kr] / mkk;
      vk = v[kr];
    }
    MPI_Bcast(&vk, 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

    if (world_rank == kw) {
      for (i = k + 1; i < size; i++) {
        m[kr * size + i] /= mkk;
      }
    }

    akj = (double *)malloc((size - k - 1) * sizeof(double));

    if (world_rank == kw) {
      for (i = k + 1; i < size; i++) {
        akj[i - k - 1] = m[kr * size + i];
      }
    }
    MPI_Bcast(akj, size - k - 1, MPI_DOUBLE, kw, MPI_COMM_WORLD);

    if (world_rank == kw) {
      for (i = 0; i < kr; i++) {
        for (j = k + 1; j < size; j++) {
          // step 5
          m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
        }
        v[i] -= m[i * size + k] * vk;
      }
      for (i = kr + 1; i < nbC; i++) {
        for (j = k + 1; j < size; j++) {
          // step 5
          m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
        }
        v[i] -= m[i * size + k] * vk;
      }
    } else {
      for (i = 0; i < nbC; i++) {
        for (j = k + 1; j < size; j++) {
          // step 5
          m[i * size + j] -= akj[j - k - 1] * m[i * size + k];
        }
        v[i] -= m[i * size + k] * vk;
      }
    }
    free(akj);
  }
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  int size;
  char *fileA, *fileB, *fileV, *fileR;
  char docstr[] = "Linear system solution with Gauss-Jordan elimination : "
                  "Ax=v\nUsage : -s size -A <path to binary file "
                  "containing A> -V <path to binary file containing v> -R "
                  "<path to binary file that will contain x>\n";
  parse_args_2mat_2vect(argc, argv, docstr, &size, &fileA, &fileB, &fileV,
                        &fileR);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // number of columns of the matrix distributed per processes
  int nbC = size / world_size;
  assert(nbC * world_size == size);

  double *v, *m;
  struct timeval ts, te, t1, t2;

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&ts, 0);
  }

  m = mat_malloc(size, world_rank, world_size);
  if (fileA == 0) {
    mat_init(m, size, world_rank, world_size);
  } else {
    mat_read_block(size, world_rank, world_size, m, fileA);
  }
  if (fileV == 0) {
    v = initVect(nbC, world_rank);
  } else {
    v = loadVect(nbC, world_rank, fileV);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  gaussJordan(size, m, v, nbC, world_size, world_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (fileR != 0) {
    saveVect(v, nbC, world_rank, fileR);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&te, 0);
    printf("%f\n",
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
    printf("%f\n",
           (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
    printf("%f\n",
           (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
  }
  if (fileB != 0) {
    mat_write_block(size, world_rank, world_size, m, fileB);
  }

  free(m);
  free(v);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
