#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct cord
{
    double value;
    struct cord *next;
};


struct vector
{
    struct vector *next;
    struct cord *cords;
};



void freeVectors(struct vector* headVec,int N)
{
    struct vector* tmp;
    int cnt = 0;
    while (cnt < N)
    {
        tmp = headVec;
        headVec = headVec->next;
        free(tmp);
        cnt += 1;
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

struct vector free_vector(struct vector *vec){
    struct vector *v_next = vec->next;
    struct cord* cordin = vec->cords;
    while (cordin != NULL){
        struct cord* next = cordin->next;
        free (cordin);
        cordin = next;
    }
    return *v_next;
}



int get_len_from_list(struct vector *vec){
    int i = 0;
    struct vector *vecs = vec;
    struct cord* cordin = vecs->cords;
    while (cordin != NULL){
        i += 1;
        cordin = cordin->next;
    }
    return i;
}


double* turn_list_to_array(struct vector *vec,int d){
    int i = 0 ;
    double * val_of_vector = (double*) calloc(d,sizeof (double*));
    struct vector *vecs = vec;
    struct cord* cordin = vecs->cords;
    while (i<d){
        val_of_vector[i]=cordin->value;
        cordin = cordin->next;
        i += 1;
    }
    return val_of_vector;
}


double **createArrayfromInput(int N ,int dim_size,struct vector head,struct cord *a){
    struct vector* v;
    double** array_of_vectors;
    double* vector_array;
    int cnt;
    array_of_vectors = (double**)calloc(N,sizeof(double*));
    cnt = 0;
    v = &head;
    while (cnt<N){
        struct vector *s;
        s = v->next;
        vector_array = turn_list_to_array(v,dim_size);
        array_of_vectors[cnt] = vector_array;
        cnt += 1;
        free_vector(v);
        v = s;
    }
    return array_of_vectors;
}



int number_of_input(struct vector * vectors){
    int i = 0;
    struct vector * ve = vectors;
    while (ve->next != NULL){
        i++;
        ve = ve->next;
    }
    return i;
}



double euclidean_distance(double *v1, double *v2, int size) {
    double dist = 0.0;
    for (int i = 0; i < size; ++i) {
        double diff = v1[i] - v2[i];
        dist += diff * diff;
    }
    return sqrt(dist);
}



double **sym(double **X, int N) {
    double **A = (double**)calloc(N, sizeof(double*));
    for (int i = 0; i < N; ++i) {
        A[i] = (double*)calloc(N, sizeof(double));
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                A[i][j] = exp( -1 * ((euclidean_distance(X[i], X[j], N)) / 2));
            }
        }
    }
    return A;
}



double **ddg(double **X, int N) {\
    double sum;
    double **A = sym(X, N);
    double **D = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; ++i) {
        D[i] = (double *)malloc(N * sizeof(double));
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



double** norm(double **X, int N){
    double **A = sym(X, N);
    double **D = ddg(X, N); // Using the ddg function
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

double **symnmf(double** X, int N){

}

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


int main(int argc, char** argv )
{
    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    int check_second_argument, N, k, iter, dim_size;
    char* first_argument;
    char*second_argument;
    double **array_of_vectors;
    double n;
    char c;
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
        printf("Invalid state's name\n");
        return 1;
    }
    int len_str = strlen(second_argument);
    if (second_argument[len_str-1]!='t' || second_argument[len_str-2]!='x' || second_argument[len_str-3]!= 't' || second_argument[len_str-4]!='.'){
        printf("Invalid file's name\n");
        return 1;
    }
    /**************************************************************************************************************/
    head_cord = malloc(sizeof(struct cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;
    head_vec = malloc(sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;
    //FILE *fptr;
    //fptr = fopen("C:\Users\yuval\CLionProjects\untitled","r");
    //if (fptr == NULL){
    // printf("Error! opening file");
    //return 1;
    //}

    while (scanf("%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
        {
            curr_cord->value = n;
            curr_vec->cords = head_cord;
            curr_vec->next = malloc(sizeof(struct vector));
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
            head_cord = malloc(sizeof(struct cord));
            curr_cord = head_cord;
            curr_cord->next = NULL;
            continue;
        }

        curr_cord->value = n;
        curr_cord->next = malloc(sizeof(struct cord));
        curr_cord = curr_cord->next;
        curr_cord->next = NULL;
    }
    /******************************************************************************************************************/
    N = number_of_input(head_vec);
    dim_size = get_len_from_list(head_vec);
    array_of_vectors = createArrayfromInput(N ,dim_size,*head_vec,head_cord);
    if(state == 1){
        create_output(sym(array_of_vectors, N),N, N);
    }
    if(state == 2){
        create_output(ddg(array_of_vectors, N),N, N);
    }
    if(state == 3){
        create_output(norm(array_of_vectors, N),N, N);
    }
    free_matrices(array_of_vectors,N);
    freeVectors(head_vec,N);
    free(head_cord);
    free(curr_vec);
    return 0;

}
