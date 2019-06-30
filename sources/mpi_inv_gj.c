#include "mpiio_dmat.h"
#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define SAVE 1

void printMat(double *mat, int size, int rank, int nprocs, char *modif) {
  int i, j;
  int nbR = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbR++;

  for (j = 0; j < nbR; j++) {
    for (i = 0; i < size; i++) {
      printf("%s%d %d %lf\n", modif, rank * nbR + j, i, mat[i + j * size]);
    }
  }
}

double *gaussJordan_inv(int size, double *m, int nprocs, int rank) {
  int i, j, k, k_div, k_loc, i_incr;
  double tmp, *kj, *B, *A;
  int nbR = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod) {
    nbR++;
    i_incr = rank * nbR;
  } else {
    i_incr = rank * nbR + mod;
  }

  B = (double *)malloc(nbR * size * sizeof(double));
  A = (double *)malloc(nbR * size * sizeof(double));

  memcpy(A, m, nbR * size * sizeof(double));

  for (i = 0; i < nbR; i++) {
    for (j = 0; j < size; j++) {
      B[i * size + j] = 0.000;
    }
    B[i * size + i + i_incr] = 1.000;
  }

  for (k = 0; k < size; k++) {
    int t = size / nprocs;
    if (k < mod * (t + 1)) {
      k_div = k / (t + 1);
      k_loc = k % (t + 1);
    } else {
      k_div = mod + (k - mod * (t + 1)) / t;
      k_loc = (k - mod * (t + 1)) % t;
    }
    // printf("r%d k%d %d %d\n", rank, k, k_div, k_loc);
    kj = (double *)malloc(2 * size * sizeof(double));

    if (rank == k_div) {
      tmp = A[k_loc * size + k];
      for (j = 0; j < size; j++) {
        B[k_loc * size + j] = B[k_loc * size + j] / tmp;
        A[k_loc * size + j] = A[k_loc * size + j] / tmp;
        kj[j] = A[k_loc * size + j];
        kj[size + j] = B[k_loc * size + j];
      }
    }

    MPI_Bcast(kj, 2 * size, MPI_DOUBLE, k_div, MPI_COMM_WORLD);

    if (rank == k_div) {
      for (i = 0; i < k_loc; i++) {
        tmp = A[i * size + k];
        for (j = 0; j < size; j++) {
          A[i * size + j] = A[i * size + j] - tmp * kj[j];
          B[i * size + j] = B[i * size + j] - tmp * kj[size + j];
        }
      }
      for (i = k_loc + 1; i < nbR; i++) {
        tmp = A[i * size + k];
        for (j = 0; j < size; j++) {
          A[i * size + j] = A[i * size + j] - tmp * kj[j];
          B[i * size + j] = B[i * size + j] - tmp * kj[size + j];
        }
      }
    } else {
      for (i = 0; i < nbR; i++) {
        tmp = A[i * size + k];
        for (j = 0; j < size; j++) {
          A[i * size + j] = A[i * size + j] - tmp * kj[j];
          B[i * size + j] = B[i * size + j] - tmp * kj[size + j];
        }
      }
    }
    free(kj);
  }

  free(A);

  return B;
}

void testLoad(int size, int world_size, int world_rank) {
  double *inv, *m;
  struct timeval ts, te, t1, t2;

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&ts, 0);
  }

  m = mat_malloc(size, world_rank, world_size);
  mat_read_block(size, world_rank, world_size, m, "a.bin");

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  inv = gaussJordan_inv(size, m, world_size, world_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  mat_write_block(size, world_rank, world_size, inv, "inv.bin");
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

  free(m);
  free(inv);
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  int size;
  if (argc == 2) {
    size = atoi(argv[1]);
  } else {
    size = 16;
  }

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  testLoad(size, world_size, world_rank);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
