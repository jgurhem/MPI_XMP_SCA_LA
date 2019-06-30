#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 16

double *initMat(int l, int rank) {
  double *mat;
  mat = (double *)malloc(l * N * sizeof(double));
  srandom((unsigned)233 * rank);
  for (int j = 0; j < l; j++) {
    for (int i = 0; i < N; i++) {
      mat[i + j * N] = 100.0 * rand() / RAND_MAX;
      ;
      // printf("a %d %d %lf\n", i, rank+l*j, mat[i*N+j]);
    }
  }
  return mat;
}

void printMat(double *mat, int nbC, int rank, char *modif) {
  for (int j = 0; j < nbC; j++) {
    for (int i = 0; i < N; i++) {
      // cyclic distribution of the columns of the matrix
      printf("%s%d %d %lf\n", modif, i, rank + j * nbC, mat[i + j * N]);
    }
  }
}

void saveMat(double *mat, int l, int rank, char *path) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype darray;
  int array_size[2] = {N, N};
  int array_distrib[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_NONE};
  int array_dargs[2] = {1, N};
  int array_psizes[2] = {world_size, 1};
  MPI_Type_create_darray(world_size /* size */, rank, 2 /* dims */, array_size,
                         array_distrib, array_dargs, array_psizes, MPI_ORDER_C,
                         MPI_DOUBLE, &darray);
  MPI_Type_commit(&darray);

  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, mat, l * N, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&darray);
}

/*
void
saveMat (double *mat, int l, int rank, char *path) {
    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
MPI_INFO_NULL, &fh); MPI_File_write_at_all (fh, rank * l * N * sizeof (double),
mat, l * N, MPI_DOUBLE, &status); MPI_File_close (&fh);
}
*/

/*
double *
initVect (int nbC, int rank) {
    double *v;
    v = (double *) malloc (nbC * sizeof (double));
    srandom ((unsigned) 21 * rank + 13456);
    for (int i = 0; i < nbC; i++) {
        v[i] = 100.0 * rand () / RAND_MAX;
        // printf("v %d %lf\n", i*nbC + rank, v[i]);
    }
    return v;
}

void
printVect (double *v, int nbC, int rank, char *modif) {
    for (int i = 0; i < nbC; i++) {
        printf ("%s%d %lf\n", modif, i * nbC + rank, v[i]);
    }
}

void
saveVect (double *v, int nbC, int rank, char *path) {
    MPI_File fh;
    MPI_Status status;
    MPI_File_open (MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
MPI_INFO_NULL, &fh); MPI_File_write_at_all (fh, rank * nbC * sizeof (double), v,
nbC, MPI_DOUBLE, &status); MPI_File_close (&fh);
}
*/

double *initVect() {
  double *v;
  v = (double *)malloc(N * sizeof(double));
  srandom((unsigned)23);
  for (int i = 0; i < N; i++) {
    v[i] = 100.0 * rand() / RAND_MAX;
    // v[i] = i;
    // printf("v %d %lf\n", i, v[i]);
  }
  return v;
}

void printVect(double *v, char *modif) {
  for (int i = 0; i < N; i++) {
    printf("%s%d %lf\n", modif, i, v[i]);
  }
}

void saveVect(double *v, char *path) {
  FILE *f;
  f = fopen(path, "w");
  for (int i = 0; i < N; i++) {
    fprintf(f, "%d %lf\n", i, v[i]);
  }
  fclose(f);
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // number of columns of the matrix distributed per processes
  int nbC = N / world_size;
  assert(nbC * world_size == N);

  double *v, *m, *aik;
  double mkk;

  if (world_rank == 0) {
    v = initVect();
    printVect(v, "v - ");
    saveVect(v, "/tmp/res/init/b1");
  }
  // saveVect (v, nbC, world_rank, "/tmp/res/init/v1");

  m = initMat(nbC, world_rank);
  // printMat (m, nbC, world_rank, "a - ");
  MPI_Barrier(MPI_COMM_WORLD);
  saveMat(m, nbC, world_rank, "/tmp/res/init/m1,1");

  MPI_Barrier(MPI_COMM_WORLD);

  int k;
  // TODO : optimize the world_size last steps were the broadcast is not
  // efficient since operations are not performed on all the processes
  for (k = 0; k < N; k++) {
    // kw : local column
    // kr : processus where the data are stored for the cyclic distribution
    int kw = k / world_size;
    int kr = k % nbC;
    if (world_rank == kr) {
      mkk = m[kw * N + k];
      // printf ("mkk %lf kr%d kw%d\n", mkk, kr, kw);
    }
    MPI_Bcast(&mkk, 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);
    /*
            if (world_rank == kr) {
                //printf ("mkk %lf vk %lf kr%d kw%d\n", mkk, v[kw], kr, kw);
                v[kw] = v[kw] / mkk;
                vk = v[kw];
            }
            MPI_Bcast (&vk, 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);
    */
    if (world_rank == 0) {
      v[k] = v[k] / mkk;
    }

    if (world_rank > kr) {
      for (int i = kw; i < nbC; i++) {
        m[i * N + k] /= mkk;
        // m[i * N + k] = -555;
      }
    } else {
      for (int i = kw + 1; i < nbC; i++) {
        m[i * N + k] /= mkk;
        // m[i * N + k] = -666;
      }
    }

    aik = (double *)malloc((N - k - 1) * sizeof(double));

    if (world_rank == kr) {
      for (int i = k + 1; i < N; i++) {
        aik[i - k - 1] = m[kw * N + i];
        // printf ("aik %lf i%d\n", aik[i - k - 1], i - k - 1);
      }
    }
    MPI_Bcast(aik, N - k - 1, MPI_DOUBLE, kr, MPI_COMM_WORLD);

    if (world_rank > kr) {
      for (int j = kw; j < nbC; j++) {
        for (int i = k + 1; i < N; i++) {
          // step 5
          m[j * N + i] -= aik[i - k - 1] * m[j * N + k];
          // m[j * N + i] = 111;
        }
        // step 4
        // v[j] -= aik[j * nbC + world_rank - k - 1] * vk;
      }
    } else {
      for (int j = kw + 1; j < nbC; j++) {
        for (int i = k + 1; i < N; i++) {
          // step 5
          m[j * N + i] -= aik[i - k - 1] * m[j * N + k];
          // m[j * N + i] = 111;
        }
        // step 4
        // v[j] -= aik[j * nbC + world_rank - k - 1] * vk;
      }
    }
    if (world_rank == 0) {
      // step 4
      for (int i = k + 1; i < N; i++) {
        v[i] -= aik[i - k - 1] * v[k];
      }
    }

    free(aik);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    printf("end\n");
  }
  // printMat (m, nbC, world_rank, "m - ");
  saveMat(m, nbC, world_rank, "/tmp/res/init/m2,2");
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Status Stat;

  for (int k = N - 1; k > 0; k--) {
    // k = 7;
    // kw : local column
    // kr : processus where the data are stored for the cyclic distribution
    int kw = k / world_size;
    int kr = k % nbC;
    // printf ("kr%d kw%d\n", kr, kw);

    if (kr == 0 && world_rank == 0) {
      for (int i = 0; i < k; i++) {
        // printf ("i%d bi %lf aik %lf  bk %lf\n", i, v[i], m[kr * N + i],
        // v[k]);
        v[i] -= m[kw * N + i] * v[k];
      }
      continue;
    }

    if (world_rank == kr) {
      aik = (double *)malloc(k * sizeof(double));
      for (int i = 0; i < k; i++) {
        aik[i] = m[kw * N + i];
      }
      MPI_Send(aik, k, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
      free(aik);
    }

    if (world_rank == 0) {
      aik = (double *)malloc(k * sizeof(double));
      MPI_Recv(aik, k, MPI_DOUBLE, kr, k, MPI_COMM_WORLD, &Stat);
      for (int i = 0; i < k; i++) {
        // printf ("i%d bi %lf aik %lf  bk %lf\n", i, v[i], aik[i], v[k]);
        v[i] -= aik[i] * v[k];
      }
      free(aik);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0) {
    printVect(v, "r - ");
    saveVect(v, "/tmp/res/init/u1");
  }
  // printVect (v, nbC, world_rank, "r - ");
  // saveVect (v, nbC, world_rank, "/tmp/res/init/v2");

  free(m);

  if (world_rank == 0)
    free(v);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
