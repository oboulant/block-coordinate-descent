#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "read_data.h"
#include "ebcd.h"

int main() {

    /*
    * Create variables
    * */
    // Read data
    double *signal;
    const char filename[255] = "./data/test2.txt";
    
    // debug
    double *res, *weights, *centered_signal;
    int i, j;

    int nb_bkps, nb_lines, nb_cols;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Set in accordance with the data file
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nb_lines = 100;
    nb_cols = 3;
    nb_bkps = 5;

    /*
    * Read data
    * */
    read_data(filename, nb_lines, nb_cols, &signal);
    printf("%f\n",signal[0]);
    printf("%f\n",signal[1]);
    printf("%f\n",signal[2]);
    printf("%f\n",signal[3]);
    printf("%f\n", signal[297]);
    printf("%f\n", signal[298]);
    printf("%f\n", signal[299]);

    res = (double*)malloc(nb_lines * nb_cols * sizeof(double));
    cumsum(signal, nb_lines, nb_cols, &res);
    printf("CCCCCCCCCCCCCCCCCCC\n");
    for (i=0 ; i<nb_lines ; i++)
    {
        for (j=0 ; j<nb_cols ; j++)
        {
            printf("%f\t", res[i*nb_cols + j]);
        }
        printf("\n");
    }
    free(res);

    res = (double*)malloc((nb_lines-1) * nb_cols * sizeof(double));
    weights = (double*)malloc((nb_lines-1) * sizeof(double));
    for (i=0 ; i<(nb_lines-1) ; i++)
    {
        weights[i] = sqrt(1.0*nb_lines/((i+1)*(nb_lines-(i+1))));
        printf("%f\t", weights[i]);
    }
    centered_signal = (double*)malloc(nb_lines * nb_cols * sizeof(double));
    center_signal(signal, nb_lines, nb_cols, &centered_signal);
    printf("centered_signal\n");
    for (i=0 ; i<nb_lines ; i++)
    {
        for (j=0 ; j<nb_cols ; j++)
        {
            printf("%f\t", centered_signal[i*nb_cols + j]);
        }
        printf("\n");
    }


    printf("leftmultiplybyXt\n");
    leftmultiplybyXt(centered_signal, nb_lines, nb_cols, weights, &res);
    for (i=0 ; i<nb_lines-1 ; i++)
    {
        for (j=0 ; j<nb_cols ; j++)
        {
            printf("%f\t", res[i*nb_cols + j]);
        }
        printf("\n");
    }
    free(res);

    res = (double*)malloc(2 * 3 * sizeof(double));
    int A_indexes[2];
    int B_indexes[3];
    A_indexes[0]= 2;
    A_indexes[1]= 4;
    B_indexes[0]= 1;
    B_indexes[1]= 5;
    B_indexes[2]= 8;
    printf("XtX\n");
    XtX(A_indexes, 2, B_indexes, 3, weights, nb_lines, &res);
    for (i=0 ; i<2 ; i++)
    {
        for (j=0 ; j<3 ; j++)
        {
            printf("%f\t", res[i*3 + j]);
        }
         printf("\n");
    }
    free(res);

    double beta[99*3];
    beta[6] = 0.5;
    beta[7] = 0.55;
    beta[8] = 0.52;
    beta[12] = 0.7;
    beta[13] = 0.77;
    beta[14] = 0.72;
    res = (double*)malloc((nb_lines-1) * nb_cols * sizeof(double));
    printf("In multiplyXtXbysparse()\n");
    multiplyXtXbysparse(A_indexes, 2, nb_lines, nb_cols, beta, weights, &res);
    for (i=0 ; i<nb_lines-1 ; i++)
    {
        for (j=0 ; j<nb_cols ; j++)
        {
            printf("%f\t", res[i*nb_cols + j]);
        }
        printf("\n");
    }
    free(res);

    res = (double*)malloc(2 * 3 * sizeof(double));
    Ebcd_Res the_res;
    ebcd(signal, nb_lines, nb_cols, 10.0, weights, -1.0, &the_res);
    printf("n_A in main is %d\n", the_res.n_A);
    for (i=0 ; i<the_res.n_A ; i++)
    {
        printf("%d\t", the_res.A[i]);
    }
    printf("\n");
    for (i=0 ; i<nb_lines ; i++)
    {
        for (j=0 ; j<nb_cols ; j++)
        {
            printf("%f\t", the_res.U[i*nb_cols + j]);
        }
        printf("\n");
    }
    free(the_res.A);
    free(the_res.U);
    fflush( stdout );
    free(res);




    





    free(signal);
    free(weights);
    free(centered_signal);


    return 0;
}