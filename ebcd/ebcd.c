#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static inline int max(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

static inline int min(int num1, int num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

void center_signal(double *signal, int n_samples, int n_dims, double **res)
{
    int i, j;
    double col_sum[n_dims];

    for (j=0 ; j<n_dims ; j++)
    {
        col_sum[j] = 0.0;
    }

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            col_sum[j] += signal[i*n_dims + j];
        }
    }

    for (j=0 ; j<n_dims ; j++)
    {
        col_sum[j] = col_sum[j] / n_samples;
    }

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = signal[i*n_dims + j] - col_sum[j];
        }
    }
}

void cumsum(double *center_signal, int n_samples, int n_dims, double **res)
{
    int i,j;

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = 0.0;
        }
    }

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] += (*res)[(i-1)*n_dims + j] + center_signal[i*n_dims + j];
        }
    }
}

void leftmultiplybyXt(double *Y, int n_samples, int n_dims, double *weights, double **res)
{
    double *Y_cumsum;

    int i, j;

    Y_cumsum = (double*)malloc(n_samples * n_dims * sizeof(double));
    cumsum(Y, n_samples, n_dims, &Y_cumsum);

    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = weights[i] * ((i+1)/n_samples*Y_cumsum[(n_samples-1)*n_dims + j] - Y_cumsum[i*n_dims + j]);
        }
    }

    free(Y_cumsum);

}

void XtX(int *A_indexes, int A_size, int *B_indexes, int B_size, double *weights, int n_samples, double **res)
{
    int i, j;
    int u, v;

    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<B_size ; j++)
        {
            u = min(A_indexes[i], B_indexes[j]);
            v = max(A_indexes[i], B_indexes[j]);
            (*res)[i*B_size + j] = weights[u] * weights[v] * (u+1.0) * (n_samples * 1.0 - (v+1.0)) / (n_samples * 1.0);
        }
    }
}

void multiplyXtXbysparse(int *A_indexes, int A_size, int n_samples, int n_dims, double *beta, double *weights, double **res)
{
    int i, j;
    double *beta_c, *s, *u, *pre_res;

    if (A_size == 0)
    {
        return;
    }

    beta_c = (double*)malloc((n_samples-1)*n_dims*sizeof(double));
    s = (double*)malloc((n_samples-1)*n_dims*sizeof(double));
    u = (double*)malloc(n_dims*sizeof(double));
    pre_res = (double*)malloc((n_samples-1)*n_dims*sizeof(double));

    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta_c[i*n_dims + j] = 0.0;
            s[i*n_dims + j] = 0.0;
            pre_res[i*n_dims + j] = 0.0;
        }
    }

    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta_c[A_indexes[i]*n_dims + j] = beta[A_indexes[i]*n_dims + j] * weights[A_indexes[i]];
        }
    }

    for (j=0 ; j<n_dims ; j++)
    {
        s[(n_samples-2)*n_dims + j] = beta_c[(n_samples-2)*n_dims + j];
    }
    for (i=1 ; i<n_samples-1 ; i++)
    {
        printf("%d\t", (n_samples-2 -i));
        printf("%f\t", beta_c[(n_samples-2 -i)*n_dims]);
        printf("%f\t", s[(n_samples-2 - (i+1))*n_dims]);
        for (j=0 ; j<n_dims ; j++)
        {   
            s[(n_samples-2 -i)*n_dims + j] = s[(n_samples-2 - i+1)*n_dims + j] + beta_c[(n_samples-2 -i)*n_dims + j];
        }
        printf("%f\n", s[(n_samples-2 -i)*n_dims]);
    }
    // Debug
    printf("s\n");
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            printf("%f\t", s[i*n_dims + j]);
        }
        printf("\n");
    }
    printf("end s\n");

    for (j=0 ; j<n_dims ; j++)
    {
        u[j] = 0.0;
    }
    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            u[j] += (A_indexes[i]*1.0+1.0)*beta_c[A_indexes[i]*n_dims + j];
        }
    }
    for (j=0 ; j<n_dims ; j++)
    {
        u[j] = u[j] / (n_samples*1.0);
    }
    // Debug
    printf("u\n");
    for (j=0 ; j<n_dims ; j++)
    {
        printf("%f\t", u[j]);
    }
    printf("\nend u\n");


    for (j=0 ; j<n_dims ; j++)
    {
        pre_res[j] = s[j] - u[j];
    }
    for (i=1 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            pre_res[i*n_dims + j] = pre_res[(i-1)*n_dims + j] + s[i*n_dims + j] - u[j];
        }
    }
    





    // for (i=0 ; i<n_samples-1 ; i++)
    // {
    //     for (j=0 ; j<n_dims ; j++)
    //     {
    //         s[i*n_dims + j] = s[i*n_dims + j] - u[j];
    //     }
    // }

    // cumsum(s, n_samples, n_dims, &pre_res);

    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = pre_res[i*n_dims + j] * weights[i];
        }
    }

    free(beta_c);
    free(s);
    free(u);
    free(pre_res);
}






// void get_C(double *d, double *center_signal, int n_samples, int n_dims, double **C)
// {
//     int i, j;
//     double *r;

//     // Allocate
//     (*C) = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
//     r = (double*)malloc(n_samples * n_dims * sizeof(double));

//     // Cumulative sum of center_signal
//     memcpy(r, center_signal, n_dims * sizeof(double));
//     for (i=1 ; i<n_samples ; i++)
//     {
//         for (j=0 ; j<n_dims ; j++)
//         {
//             r[i*n_dims + j] = r[(i-1)*n_dims + j] + center_signal[i*n_dims + j];
//         }
//     }

//     for (i=0 ; i<n_samples-1 ; i++)
//     {
//         for (j=0 ; j<n_dims ; j++)
//         {
//             (*C)[i*n_dims + j] = d[i] * (i*1.0*r[(n_samples-1)*n_dims + j] / (n_samples*1.0) - r[i*n_dims + j] );
//         }
//     }

//     free(r);
// }

// void ebcd(double *signal, int n_samples, int n_dims, double lambda)
// {

//     int i, j;
//     double *beta, *d, *X, *C;

//     /*
//      * Allocate memory
//      * */
//     beta = (double *)malloc((n_samples - 1) * n_dims * sizeof(double));

//     /*
//      * Initialize 
//      * */
//     for (i = 0 ; i < (n_samples - 1) ; i++)
//     {
//         for (j=0 ; j<n_dims ; j++)
//         {
//             beta[(i * n_dims) + j] = 0.0; 
//         }
//     }

//     get_d(n_samples, &d);
//     get_X(n_samples, d, &X);
//     center_matrix(&signal, n_samples, n_dims);
//     printf("Center Matrix\n");
//     for (i=0;i<4;i++)
//     {
//         for (j=0;j<n_dims;j++)
//         {
//             printf("%f\t", signal[i*n_dims + j]);
//         }
//         printf("\n");
//     }
//     get_C(d, signal, n_samples, n_dims, &C);
//     printf("C\n");
//     for (i=95;i<99;i++)
//     {
//         for (j=0;j<n_dims;j++)
//         {
//             printf("%f\t", C[i*n_dims + j]);
//         }
//         printf("\n");
//     }







//     free(d);
//     free(X);
//     free(C);
//     free(beta);
// }