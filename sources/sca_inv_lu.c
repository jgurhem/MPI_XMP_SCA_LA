#include "parse_args.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

void exportMat(double *mat, int size, int nb, int nc, int nr, int npcol,
               int nprow, char *path) {
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype darray;
  int array_size[2] = {size, size};
  int array_distrib[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int array_dargs[2] = {nb, nb};
  int array_psizes[2] = {nprow, npcol};
  MPI_Type_create_darray(world_size /* size */, rank, 2 /* dims */, array_size,
                         array_distrib, array_dargs, array_psizes,
                         MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
  MPI_Type_commit(&darray);

  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, mat, nc * nr, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&darray);
}

double *importMat(int size, int nb, int nc, int nr, int npcol, int nprow,
                  char *path) {
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype darray;
  int array_size[2] = {size, size};
  int array_distrib[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int array_dargs[2] = {nb, nb};
  int array_psizes[2] = {nprow, npcol};
  MPI_Type_create_darray(world_size /* size */, rank, 2 /* dims */, array_size,
                         array_distrib, array_dargs, array_psizes,
                         MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
  MPI_Type_commit(&darray);
  double *mat;
  mat = (double *)malloc(nr * nc * sizeof(double));

  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, mat, nc * nr, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&darray);
  return mat;
}

void exportVect(double *v, int size, int nb, int nrhs, int nc, int nr,
                int npcol, int nprow, char *path) {
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype darray;
  int array_size[2] = {size, nrhs};
  int array_distrib[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int array_dargs[2] = {nb, nb};
  int array_psizes[2] = {nprow, npcol};
  MPI_Type_create_darray(world_size /* size */, rank, 2 /* dims */, array_size,
                         array_distrib, array_dargs, array_psizes,
                         MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
  MPI_Type_commit(&darray);

  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, v, nc * nr, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&darray);
}

double *importVect(int size, int nb, int nrhs, int nc, int nr, int npcol,
                   int nprow, char *path) {
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_File fh;
  MPI_Status status;
  MPI_Datatype darray;
  int array_size[2] = {size, nrhs};
  int array_distrib[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
  int array_dargs[2] = {nb, nb};
  int array_psizes[2] = {nprow, npcol};
  MPI_Type_create_darray(world_size /* size */, rank, 2 /* dims */, array_size,
                         array_distrib, array_dargs, array_psizes,
                         MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
  MPI_Type_commit(&darray);

  double *v;
  v = (double *)malloc(nr * nc * sizeof(double));

  MPI_File_open(MPI_COMM_SELF, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, 0, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, v, nc * nr, MPI_DOUBLE, &status);
  MPI_File_close(&fh);
  MPI_Type_free(&darray);
  return v;
}

int main(int argc, char **argv) {

  int n;
  char *fileA, *fileI;
  char docstr[] =
      "Matrix Inversion with LU factorization\nUsage : -s size -A <path to "
      "binary file "
      "containing A> -B <path to binary file that will contain the inverse>\n";
  parse_args_2mat(argc, argv, docstr, &n, &fileA, &fileI);

  int zero = 0;
  // MPI
  int myrank, nprocs, ndims = 2, dims[2] = {0, 0};
  // MPI_Status status;
  // BLACS/SCALAPACK
  int nprow, npcol, lwork, liwork, info, ic = -1, context, myrow, mycol,
                                         desca[9];
  char order = 'R';
  // MATRIX
  int nb = 1, nr, nc, lld;
  int ia = 1, ja = 1;
  double *A, *work;
  int *ipiv, *iwork;
  struct timeval ts, te, t1, t2;

  // Initialize MPI environment
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // Initialize a default BLACS context and the processes grid
  MPI_Dims_create(nprocs, ndims, dims);
  nprow = dims[0];
  npcol = dims[1];
  blacs_get_(&ic, &zero, &context);
  blacs_gridinit_(&context, &order, &nprow, &npcol);
  blacs_gridinfo_(&context, &nprow, &npcol, &myrow, &mycol);
  // printf("w%d r%d row%d col%d\n", nprocs, myrank, myrow, mycol);

  // Computation of local matrix size
  nr = numroc_(&n, &nb, &myrow, &zero, &nprow);
  nc = numroc_(&n, &nb, &mycol, &zero, &npcol);

  // printf( "r%d nr%d nc%d ncb%d\n", myrank, nr,  nc, ncb );

  lld = max(1, nr);
  ipiv = malloc((lld + nb) * sizeof(int));
  liwork = 10;
  lwork = 10;
  work = malloc(lwork * sizeof(int));
  iwork = malloc(liwork * sizeof(int));
  lwork = liwork = -1;

  // Descriptors
  descinit_(desca, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&ts, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (fileA == NULL) {
    A = (double *)malloc(nr * nc * sizeof(double));
    srandom((unsigned)233 * myrank);
    int j;
    for (j = 0; j < nr * nc; j++) {
      A[j] = 100.0 * rand() / RAND_MAX - 50.0;
    }
  } else {
    A = importMat(n, nb, nc, nr, npcol, nprow, fileA);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // LU facto
  pdgetrf_(&n, &n, A, &ia, &ja, desca, ipiv, &info);
  pdgetri_(&n, A, &ia, &ja, desca, ipiv, work, &lwork, iwork, &liwork, &info);

  liwork = (int)iwork[0];
  lwork = (int)work[0];
  free(work);
  free(iwork);
  work = malloc(lwork * sizeof(double));
  iwork = malloc(liwork * sizeof(int));

  pdgetri_(&n, A, &ia, &ja, desca, ipiv, work, &lwork, iwork, &liwork, &info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (fileI != 0) {
    exportMat(A, n, nb, nc, nr, npcol, nprow, fileI);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&te, 0);
    printf("%f\n",
           (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0);
    printf("%f\n",
           (t2.tv_sec - ts.tv_sec) + (t2.tv_usec - ts.tv_usec) / 1000000.0);
    printf("%f\n",
           (te.tv_sec - ts.tv_sec) + (te.tv_usec - ts.tv_usec) / 1000000.0);
  }

  free(A);
  free(ipiv);

  // Close BLACS environment
  blacs_gridexit_(&context);
  blacs_exit_(&zero);

  // Close MPI environment if blacs_exit paramater is not equal zero
  // MPI_Finalize();
  return 0;
}
