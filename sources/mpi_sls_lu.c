#include "mpiio_dmat.h"
#include "parse_args.h"
#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

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

double *initVect(int size, int rank, int nprocs) {
  int nbR = size / nprocs;
  int mod = size % nprocs;
  double *v;

  if (rank < mod)
    nbR++;

  v = (double *)malloc(nbR * sizeof(double));
  srandom((unsigned)23 * rank);
  int i;
  for (i = 0; i < nbR; i++) {
    v[i] = 100.0 * rand() / RAND_MAX;
    // printf("v %d %lf\n", i, v[i]);
  }
  return v;
}

double *loadVect(int size, int rank, int nprocs, char *path) {
  int nbR = size / nprocs;
  int mod = size % nprocs;
  double *v;

  if (rank < mod)
    nbR++;

  v = (double *)malloc(nbR * sizeof(double));

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (rank < mod)
    MPI_File_read_at_all(fh, rank * nbR * sizeof(double), v, nbR, MPI_DOUBLE,
                         &status);
  else
    MPI_File_read_at_all(fh, (rank * nbR + mod) * sizeof(double), v, nbR,
                         MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  return v;
}

void saveVect(double *v, int size, int rank, int nprocs, char *path) {
  int nbR = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbR++;

  MPI_File fh;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  if (rank < mod)
    MPI_File_write_at_all(fh, rank * nbR * sizeof(double), v, nbR, MPI_DOUBLE,
                          &status);
  else
    MPI_File_write_at_all(fh, (rank * nbR + mod) * sizeof(double), v, nbR,
                          MPI_DOUBLE, &status);
  MPI_File_close(&fh);
}

void lu(int size, double *m, int nprocs, int rank) {
  int i, j, k, k_mod, k_loc, i_incr;
  double tmp, *ik;
  int nbR = size / nprocs;
  int mod = size % nprocs;

  if (rank < mod)
    nbR++;

  for (k = 0; k < size - 1; k++) {
    k_loc = k / nprocs;
    k_mod = k % nprocs;
    ik = (double *)malloc((size - k - 1) * sizeof(double));

    if (rank == k_mod) {
      for (i = k + 1; i < size; i++) {
        m[k_loc * size + i] = m[k_loc * size + i] / m[k_loc * size + k];
        ik[i - k - 1] = m[k_loc * size + i];
      }
    }

    if (rank <= k_mod) {
      k_loc++;
    }

    MPI_Bcast(ik, size - k - 1, MPI_DOUBLE, k_mod, MPI_COMM_WORLD);

    for (j = k_loc; j < nbR; j++) {
      tmp = m[j * size + k];
      for (i = k + 1; i < size; i++) {
        m[j * size + i] = m[j * size + i] - ik[i - k - 1] * tmp;
      }
    }
    free(ik);
  }
}

void resUx(int size, double *m, double *v, int nprocs, int rank) {
  int i, j, k, k_mod, k_div, js, k_loc, depla;
  double tmp, *vFull, *ik, mkk;
  int nbR = size / nprocs;
  int mod = size % nprocs;

  int counts[nprocs], displs[nprocs];
  displs[0] = 0;
  counts[nprocs - 1] = nbR;
  for (i = 0; i < nprocs - 1; i++) {
    if (i < mod) {
      counts[i] = nbR + 1;
      displs[i + 1] = displs[i] + nbR + 1;
    } else {
      counts[i] = nbR;
      displs[i + 1] = displs[i] + nbR;
    }
  }

  if (rank == 0)
    vFull = (double *)malloc(size * sizeof(double));
  if (rank < mod)
    nbR++;

  MPI_Gatherv(v, nbR, MPI_DOUBLE, vFull, counts, displs, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  for (k = size - 1; k >= 0; k--) {
    k_loc = k / nprocs;
    k_mod = k % nprocs;

    int t = size / nprocs;
    if (k < mod * (t + 1)) {
      k_div = k / (t + 1);
    } else {
      k_div = mod + (k - mod * (t + 1)) / t;
    }

    ik = (double *)malloc((k + 1) * sizeof(double));

    if (rank != 0 && k_mod != 0) {
      mkk = m[k_loc * size + k];
      MPI_Send(&m[k_loc * size], k + 1, MPI_DOUBLE, 0, k, MPI_COMM_WORLD);
    }

    if (rank == 0 && k_mod == 0) {
      vFull[k] /= m[k_loc * size + k];
      for (i = 0; i < k; i++) {
        vFull[i] -= m[k_loc * size + i] * vFull[k];
      }
    }

    if (rank == 0 && k_mod != 0) {
      MPI_Status status;
      MPI_Recv(ik, k + 1, MPI_DOUBLE, k_mod, k, MPI_COMM_WORLD, &status);
      vFull[k] /= ik[k];
      for (i = 0; i < k; i++) {
        vFull[i] -= ik[i] * vFull[k];
      }
    }
    free(ik);
  }
  MPI_Scatterv(vFull, counts, displs, MPI_DOUBLE, v, nbR, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
}

void resLx(int size, double *m, double *v, int nprocs, int rank) {
  int i, j, k, k_mod, k_div, js, k_loc, depla;
  double tmp, *jk, vk;
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

  if (rank < mod)
    nbR++;

  for (k = 0; k < size - 1; k++) {
    k_loc = k / nprocs;
    k_mod = k % nprocs;
    jk = (double *)malloc(nbR * sizeof(double));

    MPI_Scatterv(&m[k_loc * size], recvcounts, displs, MPI_DOUBLE, jk, nbR,
                 MPI_DOUBLE, k_mod, MPI_COMM_WORLD);

    int t = size / nprocs;
    if (k < mod * (t + 1)) {
      k_div = k / (t + 1);
    } else {
      k_div = mod + (k - mod * (t + 1)) / t;
    }

    if (rank == k_div)
      vk = v[k % nbR];

    MPI_Bcast(&vk, 1, MPI_DOUBLE, k_div, MPI_COMM_WORLD);

    if (k + 1 < mod * (t + 1)) {
      k_div = (k + 1) / (t + 1);
    } else {
      k_div = mod + (k + 1 - mod * (t + 1)) / t;
    }

    if (rank < k_div)
      js = nbR;
    else if (rank > k_div)
      js = 0;
    else {
      js = (k + 1) % nbR;
    }

    for (j = js; j < nbR; j++) {
      v[j] -= jk[j] * vk;
    }
    free(jk);
  }
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  int size;
  char *fileA, *fileB, *fileV, *fileR;
  char docstr[] = "Linear system solution with LU factorization: "
                  "Ax=LUx=v\nUsage : -s size -A <path to binary file "
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

  double *m, *v;
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
  if (fileV == 0) {
    v = initVect(size, world_rank, world_size);
  } else {
    v = loadVect(size, world_rank, world_size, fileV);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  lu(size, m, world_size, world_rank);
  resLx(size, m, v, world_size, world_rank);
  resUx(size, m, v, world_size, world_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (fileR != 0) {
    saveVect(v, size, world_rank, world_size, fileR);
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
    mat_write_cyclic(size, world_rank, world_size, m, fileB);
  }

  free(v);
  free(m);

  // Finalize the MPI environment.
  MPI_Finalize();
  return 0;
}
