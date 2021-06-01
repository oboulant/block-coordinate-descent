void ebcd(double *signal, int n_samples, int n_dims, double lambda);
void cumsum(double *center_signal, int n_samples, int n_dims, double **res);
void leftmultiplybyXt(double *Y, int n_samples, int n_dims, double *weights, double **res);
void center_signal(double *center_signal, int n_samples, int n_dims, double **res);
void XtX(int *A_indexes, int A_size, int *B_indexes, int B_size, double *weights, int n_samples, double **res);
void multiplyXtXbysparse(int *A_indexes, int A_size, int n_samples, int n_dims, double *beta, double *weights, double **res);
