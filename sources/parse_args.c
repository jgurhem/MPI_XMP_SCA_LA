#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void help_2mat() {
  printf("--help        print this message\n");
  printf("-s <int>      set the size of the matrix\n");
  printf("-A <string>   set file path to read A. ");
  printf("If not specified, A is generated instead\n");
  printf("-B <string>   set file path to write B.");
  printf("If not specified, B is not written\n");
}

void help_2mat_2vect() {
  printf("--help        print this message\n");
  printf("-s <int>      set the size of the matrix\n");
  printf("-A <string>   set file path to read A. ");
  printf("If not specified, A is generated instead\n");
  printf("-B <string>   set file path to write B.");
  printf("If not specified, B is not written\n");
  printf("-V <string>   set file path to read V. ");
  printf("If not specified, V is generated instead\n");
  printf("-R <string>   set file path to write R.");
  printf("If not specified, R is not written\n");
}

void parse_args_2mat(int argc, char **argv, char *docstr, int *size,
                     char **fileA, char **fileB) {
  int i;
  *size = 0;
  *fileA = 0;
  *fileB = 0;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help")) {
      printf("%s", docstr);
      help_2mat();
      continue;
    }
    if (!strcmp(argv[i], "-s")) {
      if (++i >= argc)
        continue;
      *size = atoi(argv[i]);
      continue;
    }
    if (!strcmp(argv[i], "-A")) {
      if (++i >= argc)
        continue;
      *fileA = argv[i];
      continue;
    }
    if (!strcmp(argv[i], "-B")) {
      if (++i >= argc)
        continue;
      *fileB = argv[i];
      continue;
    }
    printf("%s", docstr);
    help_2mat();
    break;
  }
}

void parse_args_2mat_2vect(int argc, char **argv, char *docstr, int *size,
                           char **fileA, char **fileB, char **fileV,
                           char **fileR) {
  int i;
  *size = 0;
  *fileA = 0;
  *fileB = 0;
  *fileV = 0;
  *fileR = 0;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help")) {
      printf("%s", docstr);
      help_2mat_2vect();
      continue;
    }
    if (!strcmp(argv[i], "-s")) {
      if (++i >= argc)
        continue;
      *size = atoi(argv[i]);
      continue;
    }
    if (!strcmp(argv[i], "-A")) {
      if (++i >= argc)
        continue;
      *fileA = argv[i];
      continue;
    }
    if (!strcmp(argv[i], "-B")) {
      if (++i >= argc)
        continue;
      *fileB = argv[i];
      continue;
    }
    if (!strcmp(argv[i], "-R")) {
      if (++i >= argc)
        continue;
      *fileR = argv[i];
      continue;
    }
    if (!strcmp(argv[i], "-V")) {
      if (++i >= argc)
        continue;
      *fileV = argv[i];
      continue;
    }
    printf("%s", docstr);
    help_2mat_2vect();
    break;
  }
}
