#include "mpiio_dmat.h"
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
  if (argc == 2) {
    n = atoi(argv[1]);
  } else {
    n = 16;
  }
  world = xmp_all_num_nodes();
  rank = xmp_node_num() - 1;
#pragma xmp template_fix t(0 : n - 1)

#pragma xmp barrier
  struct timeval ts, te, t1, t2;
#pragma xmp task on p(1)
  gettimeofday(&ts, 0);
#pragma xmp barrier
  double A0[n][n], b0[n];
  double *akj, *akj2;
  akj = malloc(n * sizeof(double));
  akj2 = malloc(n * sizeof(double));

#pragma xmp align A0[i][*] with t(i)
#pragma xmp align b0[i] with t(i)

  double btemp, bi;

  mat_read_cyclic(n, rank, world, A0, "a.bin");
#pragma xmp barrier
#pragma xmp task on p(1)
  gettimeofday(&t1, 0);
#pragma xmp barrier

#ifdef COO_OUT
  char buf[100];
  xmp_file_t *file;
  file = xmp_fopen_all("a.dat", "w");
#pragma xmp loop(i) on t(i)
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      sprintf(buf, "%d %d %lf\n", i, j, A0[i][j]);
      xmp_fwrite_shared(file, buf, strlen(buf), 1);
    }
  }
  xmp_fclose_all(file);
#pragma xmp barrier
#endif

  for (k = 0; k < n - 1; k++) {
    for (j = k; j < n; j++) {
      akj2[j] = 0.0;
    }
#pragma xmp task on t(k)
    for (j = k; j < n; j++) {
      akj2[j] = A0[k][j];
    }
    MPI_Comm comm = xmp_get_mpi_comm();
    MPI_Allreduce(akj2, akj, n, MPI_DOUBLE, MPI_SUM, comm);

#pragma xmp loop(i) on t(i)
    for (i = k + 1; i < n; i++) {
      A0[i][k] = A0[i][k] / akj[k];
      for (j = k + 1; j < n; j++) {
        A0[i][j] = A0[i][j] - A0[i][k] * akj[j];
      }
    }
  }

  free(akj);
  free(akj2);
#pragma xmp barrier
#pragma xmp task on p(1)
  gettimeofday(&t2, 0);

#pragma xmp barrier
  mat_write_cyclic(n, rank, world, A0, "lu.bin");
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
  file = xmp_fopen_all("lu.dat", "w");
#pragma xmp loop(i) on t(i)
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      sprintf(buf, "%d %d %lf\n", i, j, A0[i][j]);
      xmp_fwrite_shared(file, buf, strlen(buf), 1);
    }
  }
  xmp_fclose_all(file);
#endif

  return 0;
}
