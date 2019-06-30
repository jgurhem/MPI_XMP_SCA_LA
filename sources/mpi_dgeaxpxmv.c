#include "mpiio_dmat.h"
#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
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

void printVect(double *v, int size, int rank, int nprocs, char *modif) {
  int i;
  int nbR = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbR++;

  for (i = 0; i < nbR; i++) {
    printf("%s%d %lf\n", modif, nbR * rank + i, v[i]);
  }
}

double *loadVect(int n, int rank, int nprocs, char *path) {
  int nbR = n / nprocs;
  int mod = n % nprocs;
  double *v;

  if (rank < mod)
    nbR++;

  v = (double *)malloc(nbR * sizeof(double));

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (rank < mod) {
    MPI_File_read_at_all(fh, rank * nbR * sizeof(double), v, nbR, MPI_DOUBLE,
                         &status);
  } else {
    MPI_File_read_at_all(fh, (rank * nbR + mod) * sizeof(double), v, nbR,
                         MPI_DOUBLE, &status);
  }
  MPI_File_close(&fh);
  return v;
}

void saveVect(double *v, int n, int rank, int nprocs, char *path) {
  int nbR = n / nprocs;
  int mod = n % nprocs;

  if (rank < mod)
    nbR++;
  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  if (rank < mod)
    MPI_File_write_at_all(fh, rank * nbR * sizeof(double), v, nbR, MPI_DOUBLE,
                          &status);
  else
    MPI_File_write_at_all(fh, (rank * nbR + mod) * sizeof(double), v, nbR,
                          MPI_DOUBLE, &status);
  MPI_File_close(&fh);
}

void dgaxpxmv(int size, double *m, double *v, int nprocs, int rank) {
  int i, j, k;
  double *vFull, *vRes, tmp;
  int nbR = size / nprocs;
  int mod = size % nprocs;
  int recvcounts[nprocs], displs[nprocs];
  displs[0] = 0;
  recvcounts[nprocs - 1] = nbR;
  for (i = 0; i < nprocs - 1; i++) {
    if (i < mod) {
      recvcounts[i] = nbR + 1;
      displs[i + 1] = displs[i] + nbR + 1;
    } else {
      recvcounts[i] = nbR;
      displs[i + 1] = displs[i] + nbR;
    }
  }

  if (rank < mod) {
    nbR++;
  }

  vFull = (double *)malloc(size * sizeof(double));
  vRes = (double *)malloc(nbR * sizeof(double));
  MPI_Allgatherv(v, nbR, MPI_DOUBLE, vFull, recvcounts, displs, MPI_DOUBLE,
                 MPI_COMM_WORLD);

  for (i = 0; i < nbR; i++) {
    tmp = v[i];
    for (j = 0; j < size; j++) {
      tmp += m[i * size + j] * vFull[j];
    }
    vRes[i] = tmp;
  }

  MPI_Allgatherv(vRes, nbR, MPI_DOUBLE, vFull, recvcounts, displs, MPI_DOUBLE,
                 MPI_COMM_WORLD);

  for (i = 0; i < nbR; i++) {
    tmp = vRes[i];
    for (j = 0; j < size; j++) {
      tmp += m[i * size + j] * vFull[j];
    }
    v[i] = tmp;
  }

  free(vRes);
  free(vFull);
}

void testLoad(int size, int world_size, int world_rank) {
  double *v, *m;
  struct timeval ts, te, t1, t2;

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&ts, 0);
  }

  m = mat_malloc(size, world_rank, world_size);
  mat_read_block(size, world_rank, world_size, m, "a.bin");
  // printMat (m, size, world_rank, world_size, "a - ");

  v = loadVect(size, world_rank, world_size, "b.bin");
  // printVect (v, size, world_rank, world_size, "b - ");

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  dgaxpxmv(size, m, v, world_size, world_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  saveVect(v, size, world_rank, world_size, "r.bin");

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
  free(v);
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

  // number of columns of the matrix distributed per processes
  testLoad(size, world_size, world_rank);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
