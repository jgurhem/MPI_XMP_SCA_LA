#include "mpiio_dmat.h"
#include "parse_args.h"
#include <assert.h>
#include <fcntl.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

void printMat(double *mat, int size, int nbC, int rank, char *modif) {
  int i, j;
  for (j = 0; j < nbC; j++) {
    for (i = 0; i < size; i++) {
      // cyclic distribution of the columns of the matrix
      printf("%s%d %d %lf\n", modif, i, rank + j * nbC, mat[i + j * size]);
    }
  }
}

double *initVect(int size) {
  double *v;
  int i;
  v = (double *)malloc(size * sizeof(double));
  srandom((unsigned)23);
  for (i = 0; i < size; i++) {
    v[i] = 100.0 * rand() / RAND_MAX;
    // printf("v %d %lf\n", i, v[i]);
  }
  return v;
}

void printVect(double *v, int size, char *modif) {
  int i;
  for (i = 0; i < size; i++) {
    printf("%s%d %lf\n", modif, i, v[i]);
  }
}

double *loadVect(int size, char *path) {
  int f;
  double *v;
  v = (double *)malloc(size * sizeof(double));
  f = open(path, O_RDONLY);
  read(f, v, size * sizeof(double));
  close(f);
  return v;
}

void saveVect(double *v, int size, char *path) {
  int f;
  f = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0664);
  write(f, v, size * sizeof(double));
  close(f);
}

void gaussEliminationMPI(int size, double *m, double *v, int nbC,
                         int world_size, int world_rank) {

  MPI_Status Stat;
  int i, j, k;
  double mkk, *aik;
  // TODO : optimize the world_size last steps were the broadcast is not
  // efficient since operations are not performed on all the processes
  for (k = 0; k < size; k++) {
    // kw : local column
    // kr : processus where the data are stored for the cyclic distribution
    int kw = k / world_size;
    int kr = k % world_size;

    if (world_rank == kr) {
      mkk = m[kw * size + k];
    }
    MPI_Bcast(&mkk, 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);

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

    aik = (double *)malloc((size - k - 1) * sizeof(double));

    if (world_rank == kr) {
      for (i = k + 1; i < size; i++) {
        aik[i - k - 1] = m[kw * size + i];
      }
    }
    MPI_Bcast(aik, size - k - 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);

    if (world_rank > kr) {
      for (j = kw; j < nbC; j++) {
        for (i = k + 1; i < size; i++) {
          // step 5
          m[j * size + i] -= aik[i - k - 1] * m[j * size + k];
        }
      }
    } else {
      for (j = kw + 1; j < nbC; j++) {
        for (i = k + 1; i < size; i++) {
          // step 5
          m[j * size + i] -= aik[i - k - 1] * m[j * size + k];
        }
      }
    }
    if (world_rank == 0) {
      // step 4
      for (i = k + 1; i < size; i++) {
        v[i] -= aik[i - k - 1] * v[k];
      }
    }

    free(aik);
  }

  for (k = size - 1; k > 0; k--) {
    // k = 7;
    // kw : local column
    // kr : processus where the data are stored for the cyclic distribution
    int kw = k / world_size;
    int kr = k % world_size;

    if (kr == 0 && world_rank == 0) {
      for (i = 0; i < k; i++) {
        v[i] -= m[kw * size + i] * v[k];
      }
      continue;
    }

    if (world_rank == kr) {
      aik = (double *)malloc(k * sizeof(double));
      for (i = 0; i < k; i++) {
        aik[i] = m[kw * size + i];
      }
      MPI_Send(aik, k, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
      free(aik);
    }

    if (world_rank == 0) {
      aik = (double *)malloc(k * sizeof(double));
      MPI_Recv(aik, k, MPI_DOUBLE, kr, k, MPI_COMM_WORLD, &Stat);
      for (i = 0; i < k; i++) {
        v[i] -= aik[i] * v[k];
      }
      free(aik);
    }
  }
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  int size;
  char *fileA, *fileB, *fileV, *fileR;
  char docstr[] = "Linear system solution with Gaussian elimination : "
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
    mat_read_cyclic(size, world_rank, world_size, m, fileA);
  }
  if (world_rank == 0) {
    if (fileV == 0) {
      v = initVect(size);
    } else {
      v = loadVect(size, fileV);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  gaussEliminationMPI(size, m, v, nbC, world_size, world_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0) {
    // printVect (v,size, "r - ");
    if (fileR != 0) {
      saveVect(v, size, fileR);
    }
    gettimeofday(&te, 0);
    printf("%f\n",
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
    printf("%f\n",
           (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
    printf("%f\n",
           (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
  }

  // printMat (m, size,nbC, world_rank, "m - ");
  if (fileB != 0) {
    mat_write_cyclic(size, world_rank, world_size, m, fileB);
  }

  free(m);
  if (world_rank == 0)
    free(v);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
