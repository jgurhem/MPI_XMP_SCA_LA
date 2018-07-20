#include <string.h>
#include <stdlib.h>
#include <stdio.h>

double * inversion(double*A, int size){
    int i,j,k;
    double temp;
    double *B;

    B = malloc(size*size*sizeof(double));

    for (i=0 ; i< size; i++){
        for (j=0 ; j< size; j++){
            if(j==i)
                B[i*size+i] = 1.000;
            else
                B[i*size+j] = 0.000;
        }
    }

    for (k=0 ; k < size ; k++ ){
        temp = A[k*size+k];
        for (j=0 ; j< size ; j++){
            //printf("k%d i%d a %lf %lf\n", k, j, A[k * size + j], temp);
            A[k*size+j] = A[k*size+j]/temp;
            //printf("k%d i%d b %lf %lf\n", k, j, B[k * size + j], temp);
            B[k*size+j] = B[k*size+j]/temp;
        }

        for (i=0 ; i< size ; i++){
            if (i!=k){
                temp = A[i*size+k];
                for (j=0 ; j< size ; j++){
                    //printf("k%d i%d j%d b %lf %lf\n", k, i, j, B[i * size + j], B[k * size + j]);
                    //printf("k%d i%d j%d a %lf %lf\n", k, i, j, A[i * size + j], A[k * size + j]);
                    A[i*size+j] = A[i*size+j] - temp*A[k*size+j];
                    B[i*size+j] = B[i*size+j] - temp*B[k*size+j];
                }
            }
        }

    }
    return B;
}

double * genSingleMat(char * filePath, int bsize){
    FILE * f;
    f = fopen(filePath, "r");
    double * mat;

    mat = (double*) malloc(bsize*bsize*sizeof(double));

    int i,j,nbRow,nbCol,indI,indJ;
    double tmp;

    for(i = 0 ; i < bsize ; i++) { //row
        for(j = 0; j < bsize ; j++) { //col
            fscanf(f,"%d %d %lf",&indI,&indJ,&tmp);
            mat[indI*bsize + indJ] = tmp;
            //printf("i-%d j-%d indi-%d indj-%d %lf \n",i,j,indI,indJ,tmp);
        }
        //printf("\n");
    }
    fclose(f);

    return mat;
}

void printMatrix(double *mat, int bsize){
    int i,j;
    for(i = 0 ; i < bsize; i++) { //row
        for(j = 0; j < bsize; j++) { //col
            printf("%lf\t",mat[i*bsize + j] );
        }
        printf("\n");
    }
}

int main (int argc, char ** argv){
    if(argc != 3) {
        printf ("Wrong number of argument\n");
        exit(1);
    }
    int n = atoi (argv[2]);
    double *m;
    m = genSingleMat(argv[1], n);
    printMatrix(inversion(m, n), n);
}
