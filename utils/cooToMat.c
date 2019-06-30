#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *parseMat(char *filePath, int nbr, int nbc) {
  FILE *f;
  f = fopen(filePath, "r");
  double *mat;

  mat = (double *)malloc(nbr * nbc * sizeof(double));

  int i, j, nbRow, nbCol, indI, indJ;
  double tmp;

  for (i = 0; i < nbr; i++) {   // row
    for (j = 0; j < nbc; j++) { // col
      fscanf(f, "%d %d %lf", &indI, &indJ, &tmp);
      mat[indI * nbc + indJ] = tmp;
      // printf("i-%d j-%d indi-%d indj-%d %lf \n",i,j,indI,indJ,tmp);
    }
    // printf("\n");
  }
  fclose(f);

  return mat;
}

void printMat(double *mat, int nbr, int nbc) {
  int i, j;
  for (i = 0; i < nbr; i++) {   // row
    for (j = 0; j < nbc; j++) { // col
      printf("%lf\t", mat[i * nbc + j]);
    }
    printf("\n");
  }
}

int main(int argc, char **argv) {
  if (argc != 4)
    exit(1);

  double *m;
  int nbr = atoi(argv[2]), nbc = atoi(argv[3]);
  m = parseMat(argv[1], nbr, nbc);
  printMat(m, nbr, nbc);
  free(m);
}
