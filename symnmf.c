#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void create_output(double **vectors_array, int num_of_clusters, int dim){
    int i = 0;
    while (i < num_of_clusters){
        int j = 0;
        while (j < dim){
            double num = vectors_array[i][j];
            if (j == dim-1){
                printf("%.4f\n",num);
            }
            else {
                printf("%.4f,",num);
            }
            j += 1;
        }
        i += 1;
    }
}


void free_matrices(double ** array,int num){
    int i;
    i = 0;
    while (i <num){
        free(array[i]);
        i += 1;
    }
    free(array);
}


double euclidean_distance(double *v1, double *v2, int dim) {
    double dist = 0.0;
    for (int i = 0; i < dim; ++i) {
        double diff = v1[i] - v2[i];
        dist += diff * diff;
    }
    return dist;
}



double **sym(double **X, int N, int dim) {
    double **A = (double**)calloc(N, sizeof(double*));
    for (int i = 0; i < N; ++i) {
        A[i] = (double*)calloc(N, sizeof(double));
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                A[i][j] = exp( -1 * ((euclidean_distance(X[i], X[j], dim)) / 2));
            }
        }
    }
    return A;
}



double **ddg(double **X, int N, int dim) {\
    double sum;
    double **A = sym(X, N, dim);
    double **D = (double**)calloc(N, sizeof(double*));
    for (int i = 0; i < N; ++i) {
        D[i] = (double*)calloc(N, sizeof(double));
    }

    for (int i = 0; i < N; ++i) {
        sum = 0.0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j];
        }
        D[i][i] = sum;
    }

    // Free the memory allocated for A
    for (int i = 0; i < N; ++i) {
        free(A[i]);
    }
    free(A);

    return D;
}


double **compute_D_inv_half(double **A, int N) {
    double **D_inv_half = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; ++i) {
        D_inv_half[i] = (double *)malloc(N * sizeof(double));
        for (int j = 0; j < N; ++j) {
            D_inv_half[i][j] = (i == j) ? (1.0 / sqrt(A[i][j])) : 0.0;
        }
    }
    return D_inv_half;
}


double **matrix_multiply(double **A, double **B, int rowsA, int colsA, int colsB) {
    double **C = (double **)malloc(rowsA * sizeof(double *));
    for (int i = 0; i < rowsA; ++i) {
        C[i] = (double *)malloc(colsB * sizeof(double));
        for (int j = 0; j < colsB; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < colsA; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}




double** norm(double **X, int N, int dim){
    double **A = sym(X, N, dim);
    double **D = ddg(X, N, dim); // Using the ddg function
    double **D_inv_half = compute_D_inv_half(D, N);

    double **W = matrix_multiply(D_inv_half, A, N, N, N);
    double **W_normalized = matrix_multiply(W, D_inv_half, N, N, N);

    // Free the memory allocated for A, D, and D_inv_half
    for (int i = 0; i < N; ++i) {
        free(A[i]);
        free(D[i]);
        free(D_inv_half[i]);
    }
    free(A);
    free(D);
    free(D_inv_half);

    return W_normalized;
}

double **create_transpose(double **M,int row, int col){
    double **tr_H = (double**)calloc(col, sizeof(double*));
    for (int i = 0; i < col; ++i) {
        tr_H[i] = (double*)calloc(row, sizeof(double));}
    for (int i = 0; i<row; ++i){
        for (int j = 0; j<col;++j){
            tr_H[j][i]= M[i][j];
        }
    }
    return tr_H;
}


double **updateH(double **H ,double **W ,int iter ,double b ,double eps,int n, int k){
    int cnt = 0;
    double curr = 0;
    double dist = 2*eps;
    while ((dist >= eps) && (cnt <iter)){
        double sum = 0;
        double **WH_i = matrix_multiply(W,H,n,n,k);
        double **tr_H = create_transpose(H,n,k);
        double **H_Htr = matrix_multiply(H,tr_H,n,k,n);
        double ** H_Htr_H = matrix_multiply(H_Htr,H,n,n,k);
        for (int i = 0 ;i<n ;i++){
            for (int j =0 ; j<k ;j++){
                curr = H[i][j];
                H[i][j]=H[i][j]*((1-b) +b*((WH_i[i][j])/(H_Htr_H[i][j])));
                sum += pow((curr - H[i][j]),2);
            }
        }
        dist = sum;
        free_matrices(WH_i,n);
        free_matrices(tr_H,k);
        free_matrices(H_Htr,n);
        free_matrices(H_Htr_H,n);
        cnt += 1;

    }
    free_matrices(W,n);
    return H;
}

double **symnmf(double** X,double **H, int n,int d,int iter,double b, double eps,int k){
    double **W = norm(X,n,d);
    H = updateH(H,W,iter,b,eps,n,k);
    return H;
}


int main(int argc, char** argv )
{

    // Create a list to store the lines of the file.
    char **lines = NULL;
    int N = 0;

    // Open the file.
    FILE *fptr = fopen("C:\\Users\\yuval\\CLionProjects\\untitled\\cmake-build-debug\\input_1.txt", "r");
    if (fptr == NULL) {
        fprintf(stderr, "Error opening file.\n");
        exit(1);
    }


    // Read the lines of the file one by one.
    char line[100];
    while (fgets(line, 100, fptr) != NULL) {
        // Add the line to the list.
        lines = realloc(lines, sizeof(char *) * (N+ 1));
        lines[N] = strdup(line);
        N++;
    }

    // Close the file.
    fclose(fptr);

    // Print the list of lines.
    int dim_size = 1;
    char ch = '0';
    int j = 1;
    while (ch != 0){
        if (ch == ','){
            dim_size += 1;
        }
        ch = line[j];
        j += 1;
    }
    //creating array of vectors//
    double **array_of_vectors = (double**)calloc(N,sizeof (double*));
    for (int i = 0; i < N; ++i) {
        array_of_vectors[i] = (double*)calloc(dim_size, sizeof(double));}
    //creating array of vectors//

    for (int i = 0; i < N; i++) {
        int j = 0;
        char * token = strtok(lines[i], ",");
        while( token != NULL ) {
            array_of_vectors[i][j] = atof(token);
            token = strtok(NULL, ",");
            j += 1;
        }
    }

    // Free the memory allocated for the list.
    for (int i = 0; i < N; i++) {
        free(lines[i]);
    }
    free(lines);


    char* first_argument;
    char*second_argument;
    int state = 0;
    if (argc != 3){
        printf("An Error Has Occurred\n");
        return 1;
    }
    first_argument = argv[1];
    second_argument = argv[2];
    if (strcmp(first_argument,"sym") == 0){
        state = 1;
    }
    else if (strcmp(first_argument,"ddg") == 0){
        state = 2;
    }
    else if (strcmp(first_argument,"norm") == 0){
        state = 3;
    }
    else {
        printf("An Error Has Occurred\n");
        return 1;
    }
    int len_str = strlen(second_argument);
    if (second_argument[len_str-1]!='t' || second_argument[len_str-2]!='x' || second_argument[len_str-3]!= 't' || second_argument[len_str-4]!='.'){
        printf("An Error Has Occurred\n");
        return 1;
    }
    if(state == 1){
        create_output(sym(array_of_vectors, N,dim_size),N, N);
    }
    if(state == 2){
        create_output(ddg(array_of_vectors, N, dim_size),N, N);
    }
    if(state == 3){
        create_output(norm(array_of_vectors, N, dim_size),N, N);
    }
    
    free_matrices(array_of_vectors,N);
    return 0;

}
