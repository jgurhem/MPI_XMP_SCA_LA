#ifndef PARSE_ARGS_H
#define PARSE_ARGS_H

void help_2mat();
void help_2mat_2vect();
void parse_args_2mat(int argc, char **argv, char *docstr, int *size,
                     char **fileA, char **fileB);
void parse_args_2mat_2vect(int argc, char **argv, char *docstr, int *size,
                           char **fileA, char **fileB, char **fileV,
                           char **fileR);

#endif
