#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>   // for gettimeofday()

#include "read_data.h"
#include "ebcd.h"

int main() {

    /*
    * Create variables
    * */
    // Read data
    double *signal;
    const char filename[255] = "./data/test2.txt";
    double *weights;

    // debug
    int i;

    int nb_lines, nb_cols;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Set in accordance with the data file
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nb_lines = 100;
    nb_cols = 3;

    /*
    * Read data
    * */
    read_data(filename, nb_lines, nb_cols, &signal);

    // Weights
    weights = (double*)malloc((nb_lines-1) * sizeof(double));
    for (i=0 ; i<(nb_lines-1) ; i++)
    {
        weights[i] = sqrt(1.0*nb_lines/((i+1)*(nb_lines-(i+1))));
    }

    // run
    Ebcd_Res the_res;
    struct timeval start, end;
    gettimeofday(&start, NULL);
    ebcd_compute(signal, nb_lines, nb_cols, 10.0, weights, -1.0, &the_res);
    gettimeofday(&end, NULL);
    long seconds = (end.tv_sec - start.tv_sec);
    long micros = ((seconds * 1000000) + end.tv_usec) - (start.tv_usec); 
    printf("The elapsed time in ebcd_compute() is %ld seconds and %ld micros\n", seconds, micros);
    printf("n_A in main is %d\n", the_res.n_A);
    for (i=0 ; i<the_res.n_A ; i++)
    {
        printf("%d\t", the_res.A[i]);
    }
    printf("\n");

    free(the_res.A);
    free(the_res.U);
    free(signal);
    free(weights);


    return 0;
}