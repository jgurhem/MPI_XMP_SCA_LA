#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif
double *initMat(int size, int rank, int nprocs) {
  // block distribution of the rows
  int nbC = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbC++;

  double *mat;
  mat = (double *)malloc(nbC * size * sizeof(double));
  srandom((unsigned)233 * rank);
  int i, j;
  for (j = 0; j < nbC; j++) {
    for (i = 0; i < size; i++) {
      mat[i + j * size] = 100.0 * rand() / RAND_MAX - 50.0;
      // printf("a %d %d %lf\n", i, rank+l*j, mat[i*size+j]);
    }
  }
  return mat;
}

void saveMat(double *mat, int size, int rank, int nprocs, char *path) {
  int nbC = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbC++;

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  if (rank < mod)
    MPI_File_write_at_all(fh, rank * nbC * size * sizeof(double), mat,
                          nbC * size, MPI_DOUBLE, &status);
  else
    MPI_File_write_at_all(fh, (rank * nbC + mod) * size * sizeof(double), mat,
                          nbC * size, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
}

double *initVect(int n, int rank, int nprocs) {
  int nbC = n / nprocs;
  int mod = n % nprocs;

  if (rank < mod)
    nbC++;
  double *v;
  int i;
  v = (double *)malloc(nbC * sizeof(double));
  srandom((unsigned)21 * rank + 13456);
  for (i = 0; i < nbC; i++) {
    v[i] = 100.0 * rand() / RAND_MAX - 50.0;
    // printf("v %d %lf\n", i*nbC + rank, v[i]);
  }
  return v;
}

void saveVect(double *v, int n, int rank, int nprocs, char *path) {
  int nbC = n / nprocs;
  int mod = n % nprocs;

  if (rank < mod)
    nbC++;
  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  if (rank < mod)
    MPI_File_write_at_all(fh, rank * nbC * sizeof(double), v, nbC, MPI_DOUBLE,
                          &status);
  else
    MPI_File_write_at_all(fh, (rank * nbC + mod) * sizeof(double), v, nbC,
                          MPI_DOUBLE, &status);
  MPI_File_close(&fh);
}

int main(int argc, char **argv) {

  int n;
  if (argc == 2) {
    n = atoi(argv[1]);
  } else {
    n = 16;
  }

  MPI_Init(NULL, NULL);

  int nprocs, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  double *v, *m;

  v = initVect(n, world_rank, nprocs);
  m = initMat(n, world_rank, nprocs);

  saveVect(v, n, world_rank, nprocs, "b.bin");
  saveMat(m, n, world_rank, nprocs, "a.bin");

  free(v);
  free(m);

  MPI_Finalize();
  return 0;
}
