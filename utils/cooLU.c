#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void factLU(double *A, int size) {

  int i, j, k;
  for (k = 0; k < size - 1; k++) {
    for (i = k + 1; i < size; i++) {
      A[i * size + k] = A[i * size + k] / A[k * size + k];
    }
    for (i = k + 1; i < size; i++) {
      for (j = k + 1; j < size; j++) {
        A[i * size + j] = A[i * size + j] - A[i * size + k] * A[k * size + j];
      }
    }
  }
}

double *genSingleMat(char *filePath, int bsize) {
  FILE *f;
  f = fopen(filePath, "r");
  double *mat;

  mat = (double *)malloc(bsize * bsize * sizeof(double));

  int i, j, nbRow, nbCol, indI, indJ;
  double tmp;

  for (i = 0; i < bsize; i++) {   // row
    for (j = 0; j < bsize; j++) { // col
      fscanf(f, "%d %d %lf", &indI, &indJ, &tmp);
      mat[indI * bsize + indJ] = tmp;
      // printf("i-%d j-%d indi-%d indj-%d %lf \n",i,j,indI,indJ,tmp);
    }
    // printf("\n");
  }
  fclose(f);

  return mat;
}

void printU(double *mat, int bsize) {
  printf("u\n");
  int i, j;
  for (i = 0; i < bsize; i++) {   // row
    for (j = 0; j < bsize; j++) { // col
      if (i <= j)
        printf("%lf\t", mat[i * bsize + j]);
      else
        printf("%lf\t", 0.0);
    }
    printf("\n");
  }
}

void printL(double *mat, int bsize) {
  printf("l\n");
  int i, j;
  for (i = 0; i < bsize; i++) {   // row
    for (j = 0; j < bsize; j++) { // col
      if (i > j)
        printf("%lf\t", mat[i * bsize + j]);
      else if (i == j)
        printf("%lf\t", 1.0);
      else
        printf("%lf\t", 0.0);
    }
    printf("\n");
  }
}

void printMatrix(double *mat, int bsize) {
  int i, j;
  for (i = 0; i < bsize; i++) {   // row
    for (j = 0; j < bsize; j++) { // col
      printf("%lf\t", mat[i * bsize + j]);
    }
    printf("\n");
  }
}

int main(int argc, char **argv) {
  if (argc != 3) {
    printf("Wrong number of argument\n");
    exit(1);
  }
  int n = atoi(argv[2]);
  double *m;
  m = genSingleMat(argv[1], n);
  factLU(m, n);
  printMatrix(m, n);
  // printU(m, n);
  // printL(m, n);
  free(m);
}
