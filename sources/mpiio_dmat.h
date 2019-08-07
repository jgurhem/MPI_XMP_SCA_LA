#ifndef MPIIO_DMAT_H
#define MPIIO_DMAT_H

double *mat_malloc(int size, int rank, int nprocs);
void mat_init(double *mat, int size, int rank, int nprocs);
void mat_read_cyclic(int size, int rank, int nprocs, void *mat, char *path);
void mat_write_cyclic(int size, int rank, int nprocs, void *mat, char *path);

double *vect_malloc(int size, int rank, int nprocs);
void vect_init(double *v, int size, int rank, int nprocs);
void vect_read_cyclic(int size, int rank, int nprocs, void *vect, char *path);
void vect_write_cyclic(int size, int rank, int nprocs, void *vect, char *path);

void mat_read_block(int size, int rank, int nprocs, void *mat, char *path);
void mat_write_block(int size, int rank, int nprocs, void *mat, char *path);

#endif
