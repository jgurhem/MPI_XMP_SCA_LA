#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "wrong number of arguments\n");
    exit(1);
  }

  FILE *in, *out;
  in = fopen(argv[1], "r");
  if (in == NULL) {
    printf("open error : %s\n", argv[1]);
    exit(1);
  }
  out = fopen(argv[2], "w");
  if (out == NULL) {
    printf("open error : %s\n", argv[2]);
    exit(1);
  }
  double d;
  int nbRow = atoi(argv[3]), nbCol = atoi(argv[4]);
  for (int j = 0; j < nbCol; j++) {
    for (int i = 0; i < nbRow; i++) {
      int rt;
      rt = fread(&d, sizeof(double), 1, in);
      if (rt > 0) {
        fprintf(out, "%d %d %lf\n", i, j, d);
        // fprintf (stdout, "%d %d %lf\n", i, j ,d);
      } else {
        exit(0);
      }
    }
  }
}
