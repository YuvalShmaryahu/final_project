# ifndef SYMNMF_H_
# define SYMNMF_H_

double **symnmf(double** X,double **H, int n,int d,int iter,double b, double eps,int k);
double** norm(double **X, int N, int dim);
double **ddg(double **X, int N, int dim);
double **sym(double **X, int N, int dim);
void create_output(double **vectors_array, int num_of_clusters, int dim);

#endif
