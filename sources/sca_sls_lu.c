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

  int n, nrhs = 1;
  if (argc == 2) {
    n = atoi(argv[1]);
  } else {
    n = 16;
  }

  int zero = 0, uno = 1;
  // MPI
  int myrank, nprocs, ndims = 2, dims[2] = {0, 0};
  // MPI_Status status;
  // BLACS/SCALAPACK
  int nprow, npcol, info, ic = -1, context, myrow, mycol, desca[9], descb[9];
  char order = 'R';
  // MATRIX
  int nb = 1, nr, nc, ncb, lld;
  double *A, *b;
  int *ipiv;
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
  ncb = numroc_(&nrhs, &nb, &mycol, &zero, &npcol);

  // printf( "r%d nr%d nc%d ncb%d\n", myrank, nr,  nc, ncb );

  lld = max(1, nr);
  ipiv = malloc((lld + nb) * sizeof(int));

  // Descriptors
  descinit_(descb, &n, &nrhs, &nb, &nb, &zero, &zero, &context, &lld, &info);
  descinit_(desca, &n, &n, &nb, &nb, &zero, &zero, &context, &lld, &info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&ts, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  A = importMat(n, nb, nc, nr, npcol, nprow, "a.bin");
  b = importVect(n, nrhs, nb, ncb, nr, npcol, nprow, "b.bin");

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&t1, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Linear system equations solver
  pdgesv_(&n, &uno, A, &uno, &uno, desca, ipiv, b, &uno, &uno, descb, &info);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    gettimeofday(&t2, 0);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  exportVect(b, n, nrhs, nb, ncb, nr, npcol, nprow, "r.bin");

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
  free(b);
  free(ipiv);

  // Close BLACS environment
  blacs_gridexit_(&context);
  blacs_exit_(&zero);

  // Close MPI environment if blacs_exit paramater is not equal zero
  // MPI_Finalize();
  return 0;
}
