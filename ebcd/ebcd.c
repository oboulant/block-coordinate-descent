#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEFAULT_GAIN_TOL 1e-8
#define DEFAULT_BETA_TOL 1e-12
//#define DEFAULT_MAX_IT_OUTER 1e5
#define DEFAULT_MAX_IT_OUTER 100
//#define DEFAULT_MAX_IT 1e5
#define DEFAULT_MAX_IT 100


static inline int max_i(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

static inline int min_i(int num1, int num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

static inline double max_f(double num1, double num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

static inline double min_f(double num1, double num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

static inline void remove_element_i(int *array_from, int n_el, int i_remove, int **array_to)
{
    int i, c;
    c = 0;

    for (i=0 ; i<n_el ; i++)
    {
        if (i == i_remove)
            continue;
        (*array_to)[c] = array_from[i];
        c++;
    }
}

static inline double frobenius_norm(double *Y, int n_samples, int n_dims)
{
    int i, j;
    double res;

    res = 0.0;

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            res += pow(Y[i*n_dims + j], 2);
        }
    }

    return sqrt(res);
}

static inline int max_index_array(double *Y, int n)
{
    int i, res;
    double max;

    if (n<1)
    {
        return -1;
    }

    res = 0;
    max = Y[0];

    if (n == 1)
    {
        return res;
    }

    for (i=1 ; i<n ; i++)
    {
        if (Y[i] > max)
        {
            res = i;
            max = Y[i];
        }
    }
    return res;
}

static inline int check_i_in_A_remove_inplace(int **A, int array_size, int val)
{
    int i, j, k;

    for (i=0 ; i<array_size ; i++)
    {
        if ((*A)[i] == val)
        {
            for (j=i ; j<array_size-1 ; j++)
            {
                (*A)[j] = (*A)[j+1];
            }
            (*A)[array_size-1] = '\0';
            return 1;
        }
    }
    return 0;
}

static inline void add_element_in_A_inplace(int **A, int max_array_size, int array_size, int val)
{

    (*A)[array_size] = val;
    if (array_size == max_array_size)
    {
        return;
    }
    (*A)[array_size+1] = '\0';
    return;
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
        col_sum[j] = col_sum[j] / (n_samples*1.0);
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
            u = min_i(A_indexes[i], B_indexes[j]);
            v = max_i(A_indexes[i], B_indexes[j]);
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
        for (i=0 ; i<n_samples-1 ; i++)
        {
            for (j=0 ; j<n_dims ; j++)
            {
                (*res)[i*n_dims + j] = 0.0;
            }
        }

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
        // printf("%d\t", (n_samples-2 -i));
        // printf("%f\t", beta_c[(n_samples-2 -i)*n_dims]);
        // printf("%f\t", s[(n_samples-2 - (i+1))*n_dims]);
        for (j=0 ; j<n_dims ; j++)
        {   
            s[(n_samples-2 -i)*n_dims + j] = s[(n_samples-2 - i+1)*n_dims + j] + beta_c[(n_samples-2 -i)*n_dims + j];
        }
        // printf("%f\n", s[(n_samples-2 -i)*n_dims]);
    }
    // Debug
    // printf("s\n");
    // for (i=0 ; i<n_samples-1 ; i++)
    // {
    //     for (j=0 ; j<n_dims ; j++)
    //     {
    //         printf("%f\t", s[i*n_dims + j]);
    //     }
    //     printf("\n");
    // }
    // printf("end s\n");

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
    // printf("u\n");
    // for (j=0 ; j<n_dims ; j++)
    // {
    //     printf("%f\t", u[j]);
    // }
    // printf("\nend u\n");


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


void ebcd(double *signal, int n_samples, int n_dims, double lambda, double *weights, double tol, double **res)
{
    int i, j, n_A, p, q, n_A_not_indexes;
    int it, A_idx;
    int *A, *Ai, *A_not_indexes;
    int global_sol, max_it, it_counter, i_max_norm_A_not_indexes;
    double tol_c, lagr;
    double XitX_dot_beta_Ai, gammai;
    double *centered_signal, *beta, *C;
    double *gain, *S_i, *XitX, *new_beta_i, *S, *normS;
    double temp_d, max_norm_A_not_indexes;
    double *temp_d_array;

    tol_c = tol;
    if (tol < 0.0)
    {
        tol_c = DEFAULT_GAIN_TOL;
    }

    /*
    *   Center signal
    */
    centered_signal = (double*)malloc(n_samples * n_dims * sizeof(double));
    center_signal(signal, n_samples, n_dims, &centered_signal);

    /*
    * Initialize A and beta
    */
    n_A = 0;
    A = (int*)malloc((n_samples-1) * sizeof(int));
    A_not_indexes = (int*)malloc((n_samples-1) * sizeof(int));
    A[0] = '\0';
    beta = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta[i*n_dims + j] = 0.0;
        }
    }

    /*
    * Compute C = X'*Y
    */
    C = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
    leftmultiplybyXt(centered_signal, n_samples, n_dims, weights, &C);
    // Debug
    printf("##################### C\n");
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            printf("%f\t", C[i*n_dims + j]);
        }
        printf("\n");
    }

    global_sol = 0;

    it_counter = 0;

    while (global_sol == 0 && it_counter < DEFAULT_MAX_IT_OUTER)
    {
        // Debug
        printf("#### Iteration %d\n", it_counter);
        it_counter+=1;

        gain = (double*)malloc(n_A * sizeof(double));
        for (i=0 ; i<n_A ; i++)
        {
            gain[i] = 2.0*tol_c;
        }
        // Debug
        printf("#### gain\n");
        for (i=0 ; i<n_A ; i++)
        {
            printf("%.10f\t", gain[i]);
        }
        printf("\n");
        printf("#### beta\n");
        for (i=0 ; i<n_A ; i++)
        {
            for (p=0 ; p<n_dims ; p++)
            {
                printf("%f\t", beta[A[i]*n_dims + p]);
            }
             printf("\n");
        }

        Ai = (int*)malloc((n_A-1) * sizeof(int));
        S_i = (double*)malloc(n_dims * sizeof(double));
        new_beta_i = (double*)malloc(n_dims * sizeof(double));
        temp_d_array = (double*)malloc(1 * n_dims * sizeof(double));
        XitX = (double*)malloc(1 * (n_A-1) * sizeof(double));
        for (it=0 ; it<DEFAULT_MAX_IT ; it++)
        {
            // debug
            printf("###### it %d\n", it);
            if (n_A == 0)
            {
                printf("###### Breaking beta optim on n_A == 0\n");
                break;
            }
            // Update beta
            A_idx = it % n_A;
            i = A[A_idx];
            remove_element_i(A, n_A, A_idx, &Ai);
            // Debug
            printf("A_idx : %d\n", A_idx);
            printf("i : %d\n", i);
            // debug
            printf("###### A with n_A=%d\n", n_A);
            for (p=0 ; p<n_A ; p++)
            {
                printf("%d\t", A[p]);   
            }
            printf("\n");
            printf("######Ai\n");
            for (p=0 ; p<n_A-1 ; p++)
            {
                printf("%d\t", Ai[p]);   
            }
            printf("\n");
            
            for (p=0 ; p<n_dims ; p++)
            {
                S_i[p] = C[i*n_dims + p];
            }
            // Debug
            printf("S_i before\n");
            for (p=0 ; p<n_dims ; p++)
            {
                printf("%f\t", S_i[p]); 
            }
            printf("\n");
            // Reprendre ici !!!!
            if (n_A > 1)
            {
                XtX(&i, 1, Ai, n_A-1, weights, n_samples, &XitX);
                printf("XitX\n");
                for (q=0 ; q<n_A-1 ; q++)
                {
                    printf("%f\t", XitX[q]);
                }
                printf("\n");
                for (p=0 ; p<n_dims ; p++)
                {
                    XitX_dot_beta_Ai = 0.0;
                    for (q=0 ; q<n_A-1 ; q++)
                    {
                        XitX_dot_beta_Ai += XitX[q] * beta[Ai[q]*n_dims + p];
                    }
                    S_i[p] -= XitX_dot_beta_Ai;
                }
            }

            // Debug
            printf("S_i\n");
            for (p=0 ; p<n_dims ; p++)
            {
                printf("%f\t", S_i[p]); 
            }
            printf("\n");
            
            gammai = (i + 1.0) * (n_samples - i - 1.0) * pow(weights[i], 2) / (n_samples*1.0);
            // Debug
            printf("gammai : %f\n", gammai);
            temp_d = max_f(1.0 - lambda / frobenius_norm(S_i, 1, n_dims), 0.0);
            // Debug
            printf("temp_d : %f\n", temp_d);
            // Debug
            printf("frobenius_norm(S_i, 1, n_dims) : %f\n", frobenius_norm(S_i, 1, n_dims));
            for (p=0 ; p<n_dims ; p++)
            {
                new_beta_i[p] = temp_d * S_i[p] / gammai;
                temp_d_array[p] = beta[i*n_dims + p] - new_beta_i[p];
                // debug
                printf("new_beta_i[p]  %f\t", new_beta_i[p] );
                printf("temp_d_array[%d] : %f\t", p, temp_d_array[p]);
            }
             // debug
            printf("\n");
            gain[A_idx] = frobenius_norm(temp_d_array, 1, n_dims);
            for (p=0 ; p<n_dims ; p++)
            {
                beta[i*n_dims + p] = new_beta_i[p];
            }
            // Debug
            printf("gain\n");
            for (p=0 ; p<n_A ; p++)
            {
                printf("%.12f\t", gain[p]);
            }
            printf("\n");
            printf("betai\n");
            for (p=0 ; p<n_dims ; p++)
            {
                printf("%f\t", beta[i*n_dims + p]);
            }
            printf("\n");

            printf("max_index_array(gain, n_A) is %d\n", max_index_array(gain, n_A));
            printf("gain[max_index_array(gain, n_A)] is %f\n", gain[max_index_array(gain, n_A)]);

            if (gain[max_index_array(gain, n_A)] < tol_c)
            {
                // Debug
                printf("In breaking\n");
                break;
            }
        }
        // Check where it should be called, just after the for loop, or before going back to the while loop
        printf("bla1\n");
        free(Ai);
        printf("bla2\n");
        free(S_i);
        printf("bla3\n");
        free(new_beta_i);
        printf("bla4\n");
        free(temp_d_array);
        printf("bla5\n");
        free(XitX);
        printf("bla6\n");
        free(gain);
        printf("bla7\n");

        // Remove from active set the zero coefficients
        for (i=(n_samples-2) ; i>=0 ; i--)
        {
            temp_d = 0.0;
            for (j=0 ; j<n_dims ; j++)
            {
                temp_d += fabs(beta[i*n_dims + j]);
            }
            if (temp_d < DEFAULT_BETA_TOL)
            {
                for (j=0 ; j<n_dims ; j++)
                {
                    beta[i*n_dims + j] = 0.0;
                }
                if (check_i_in_A_remove_inplace(&A, n_A, i) == 1)
                {
                    n_A -= 1;
                }
            }
        }
        printf("After removing from active set : n_A = %d\n", n_A);

        // Check optimality
        temp_d_array = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
        S = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
        normS = (double*)malloc((n_samples-1) * sizeof(double));
        // Debug
        printf("Debug\n");
        printf("n_A is %d\n", n_A);
        multiplyXtXbysparse(A, n_A, n_samples, n_dims, beta, weights, &temp_d_array);
        // print S
        printf("S\n");
        for (i=0 ; i<n_samples-1 ; i++)
        {
            normS[i] = 0.0;
            for (j=0 ; j<n_dims ; j++)
            {
                S[i*n_dims + j] = C[i*n_dims + j] - temp_d_array[i*n_dims + j];
                printf("%f\t", S[i*n_dims + j]);
                normS[i] += pow(S[i*n_dims + j ], 2);
            }
            printf("\n");
        }
        // Debug
        printf("normS\n");
        for (i=0 ; i<n_samples-1 ; i++)
        {
            printf("%f\n", normS[i]);
        }
        printf("bla8\n");
        free(temp_d_array);
        printf("bla9\n");
        free(S);
        printf("bla10\n");

        if (n_A > 0)
        {
            // Debug
            printf("Check optimality with n_A > 0\n");
            // At optimality we must have normS(i)=lambda^2 for i in AS and
            // normS(i)<lambda^2 for i not in AS.
            lagr =0.0;
            for (i=0 ; i<n_A ; i++)
            {
                if (normS[A[i]] > lagr)
                {
                    lagr = normS[A[i]];
                }
            }
            // Debug
            printf("lagr is %f\n", lagr);
            if (pow(lambda, 2) < lagr)
            {
                lagr = pow(lambda, 2);
            }
            // Debug
            printf("pow(lambda, 2) is %f\n", pow(lambda, 2));
            printf("lagr is %f\n", lagr);

            n_A_not_indexes = n_samples-1;
            for (i=0 ; i<n_samples-1 ; i++)
            {
                A_not_indexes[i] = i;
            }
            for (i=0 ; i<n_A ; i++)
            {
                if(check_i_in_A_remove_inplace(&A_not_indexes, n_A_not_indexes, A[i]) == 1)
                {
                    // Should be the case
                    // Debug
                    printf("######################################## It is the case\n");
                    n_A_not_indexes -= 1;
                }
            }
            // Debug
            printf("n_A_not_indexes is %d\n", n_A_not_indexes);
            for (i=0 ; i<n_A_not_indexes ; i++)
            {
                printf("%d\t", A_not_indexes[i]);
            }
            printf("\n");            
            max_norm_A_not_indexes = 0.0;
            i_max_norm_A_not_indexes = -1;
            for (i=0 ; i<n_A_not_indexes ; i++)
            {
                if (normS[A_not_indexes[i]] > max_norm_A_not_indexes)
                {
                    max_norm_A_not_indexes = normS[A_not_indexes[i]];
                    i_max_norm_A_not_indexes = A_not_indexes[i];
                }
            }
            // Debug
            printf("max_norm_A_not_indexes is %f\n", max_norm_A_not_indexes);
            printf("i_max_norm_A_not_indexes is %d\n", i_max_norm_A_not_indexes);
            if ((n_A_not_indexes == 0) || (max_norm_A_not_indexes < lagr + tol_c))
            {
                // Optimality conditions are fulfilled: we have found the global
                // solution
                global_sol = 1;
            }
            else
            {
                // Otherwise we add the block that violates most the optimality
                // condition
                add_element_in_A_inplace(&A, n_samples-1, n_A, i_max_norm_A_not_indexes);
                n_A += 1;
            }
        }
        else
        {
            // Debug
            printf("#### In case n_A == 0\n");
            i = max_index_array(normS, n_samples);
            if (normS[i] < pow(lambda, 2) + tol_c)
            {
                printf("###### Setting global_sol to 1\n");
                global_sol = 1;
            }
            else
            {
                printf("###### Adding %d in A\n", i);
                add_element_in_A_inplace(&A, n_samples-1, n_A, i);
                n_A = 1;
                for (j=0 ; j<n_dims ; j++)
                {
                    beta[i*n_dims + j] = 0.0;
                }
            }
        }
        printf("bla11\n");
        free(normS);
        printf("bla12\n");
        // if (it_counter == 5)
        // {
        //     return;
        // }
    }
    
    // Debug
    printf("##################### beta\n");
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            printf("%f\t", beta[i*n_dims + j]);
        }
        printf("\n");
    }
    printf("##################### A\n");
    for (i=0 ; i<n_A ; i++)
    {
        printf("%d\t", A[i]);
    }
    printf("\n");

    // Reconstruct U
    // TODO
    printf("bla13\n");
    free(centered_signal); 
    printf("bla14\n");
    free(C);
    printf("bla15\n");
    free(beta);
    printf("bla16\n");
    free(A_not_indexes);
    printf("bla17\n");
    free(A);
    // FREE A ?
    return;
}
