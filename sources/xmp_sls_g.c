#include "mpiio_dmat.h"
#include "parse_args.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <xmp.h>

#ifdef COO_OUT
#include <string.h>
#include <xmp_io.h>
#endif

#pragma xmp nodes p(*)
#pragma xmp template t( :)
#pragma xmp distribute t(cyclic) onto p

int main(int argc, char **argv) {
  int i, j, k, n, world, rank;
  char *fileA0, *fileB0, *fileV, *fileR;
  parse_args_2mat_2vect(argc, argv, &n, &fileA0, &fileB0, &fileV, &fileR);
  world = xmp_all_num_nodes();
  rank = xmp_node_num() - 1;
#pragma xmp template_fix t(0 : n - 1)

  struct timeval ts, te, t1, t2;
#pragma xmp barrier
#pragma xmp task on p(1)
  gettimeofday(&ts, 0);
#pragma xmp barrier
  double A0[n][n], b0[n];
  double *akj, *akj2;
  akj = malloc(n * sizeof(double));
  akj2 = malloc(n * sizeof(double));

#pragma xmp align A0[i][*] with t(i)
#pragma xmp align b0[i] with t(i)

  double btemp, akk;

  if (fileA0 == 0) {
    mat_init(A0, n, rank, world);
  } else {
    mat_read_cyclic(n, rank, world, A0, fileA0);
  }
  if (fileV == 0) {
    vect_init(b0, n, rank, world);
  } else {
    vect_read_cyclic(n, rank, world, b0, fileV);
  }
#pragma xmp barrier

#ifdef COO_OUT
  char buf[100];
  xmp_file_t *fileA, *fileB;
  fileA = xmp_fopen_all("a.dat", "w");
  fileB = xmp_fopen_all("v.dat", "w");
#pragma xmp loop(i) on t(i)
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      sprintf(buf, "%d %d %lf\n", i, j, A0[i][j]);
      xmp_fwrite_shared(fileA, buf, strlen(buf), 1);
    }
    sprintf(buf, "%d %d %lf\n", i, 0, b0[i]);
    xmp_fwrite_shared(fileB, buf, strlen(buf), 1);
  }
  xmp_fclose_all(fileA);
  xmp_fclose_all(fileB);
#endif

#pragma xmp task on p(1)
  gettimeofday(&t1, 0);
#pragma xmp barrier

  for (k = 0; k < n; k++) {

    akk = 0;
    btemp = 0;
#pragma xmp task on t(k)
    {
      akk = A0[k][k];
      btemp = b0[k];
    }
#pragma xmp reduction(+ : akk, btemp)

#pragma xmp loop(i) on t(i)
    for (i = k + 1; i < n; i++) {
      A0[i][k] = A0[i][k] / akk;
      b0[i] = b0[i] - A0[i][k] * btemp;
    }

    for (j = k + 1; j < n; j++) {
      akj2[j] = 0.0;
    }

#pragma xmp task on t(k)
    for (j = k + 1; j < n; j++) {
      akj2[j] = A0[k][j];
    }
    MPI_Comm comm = xmp_get_mpi_comm();
    MPI_Allreduce(akj2, akj, n, MPI_DOUBLE, MPI_SUM, comm);

#pragma xmp loop(i) on t(i)
    for (i = k + 1; i < n; i++) {
      for (j = k + 1; j < n; j++) {
        A0[i][j] = A0[i][j] - A0[i][k] * akj[j];
      }
    }
  }

  for (k = n - 1; k >= 0; k--) {
    btemp = 0;
#pragma xmp task on t(k)
    {
      b0[k] = b0[k] / A0[k][k];
      btemp = b0[k];
    }
#pragma xmp reduction(+ : btemp)

#pragma xmp loop(i) on t(i)
    for (i = 0; i < k; i++) {
      b0[i] = b0[i] - A0[i][k] * btemp;
    }
  }
  free(akj);
  free(akj2);
#pragma xmp barrier
#pragma xmp task on p(1)
  gettimeofday(&t2, 0);
#pragma xmp barrier

  if (fileR != 0)
    vect_write_cyclic(n, rank, world, b0, fileR);
#pragma xmp barrier
#pragma xmp task on p(1)
  {
    gettimeofday(&te, 0);
    printf("%f\n",
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
    printf("%f\n",
           (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
    printf("%f\n",
           (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
  }

#ifdef COO_OUT
  fileB = xmp_fopen_all("r.dat", "w");
#pragma xmp loop(i) on t(i)
  for (i = 0; i < n; i++) {
    sprintf(buf, "%d %d %lf\n", i, 0, b0[i]);
    xmp_fwrite_shared(fileB, buf, strlen(buf), 1);
  }
  xmp_fclose_all(fileB);
#endif

  return 0;
}
