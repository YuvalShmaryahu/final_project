# ifndef SYMNMF_H_
# define SYMNMF_H_

double **symnmf(double** X,double **H, int n,int d,int iter,double b, double eps,int k);
double** norm(double **X, int N, int dim);
double **ddg(double **X, int N, int dim);
double **sym(double **X, int N, int dim);

#endif
