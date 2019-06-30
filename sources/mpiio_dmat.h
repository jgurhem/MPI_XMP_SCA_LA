#ifndef MPIIO_DMAT_H
#define MPIIO_DMAT_H

double *mat_malloc(int size, int rank, int nprocs);
void mat_read_cyclic(int size, int rank, int nprocs, void *mat, char *path);
void mat_write_cyclic(int size, int rank, int nprocs, void *mat, char *path);

void vect_read_cyclic(int size, int rank, int nprocs, void *vect, char *path);
void vect_write_cyclic(int size, int rank, int nprocs, void *vect, char *path);

void mat_read_block(int size, int rank, int nprocs, void *mat, char *path);
void mat_write_block(int size, int rank, int nprocs, void *mat, char *path);

#endif
