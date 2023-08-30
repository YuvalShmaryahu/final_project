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
                sum += (curr - H[i][j])*(curr - H[i][j]);
            }
        }
        dist = sqrt(sum);
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
    int num_of_clusters =3;
    /**********Initalizing H***********/
    double F[10][3] = {{0.1915330298744675, 0.24959733187150154, 0.210361251846407}, {0.19016136851628856, 0.14785329944976214, 0.22541365248570946}, {0.15271563802220298, 0.3112240926612766, 0.33631324115595573}, {0.13381907579681968, 0.2763078792454998, 0.18458154863277862}, {0.19824456780210784, 0.3230283644720431, 0.02479121114844506}, {0.03040766789311094, 0.007056114496527709, 0.29057995222458244}, {0.271572618312309, 0.30362966943935477, 0.3415326606165121}, {0.2789021408169059, 0.16105387320889544, 0.2724005822920931}, {0.041277153315709414, 0.2233290748462612, 0.05002954425034687}, {0.32968448956878926, 0.18212232294658148, 0.14471483877795904}};
    double **H = (double**)calloc(10, sizeof(double*));
    for (int i =0;i<10;i++){
        H[i]=F[i];
    }
    /**********Initalizing H***********/
    array_of_vectors = createArrayfromInput(N ,dim_size,*head_vec,head_cord);
    double eps  = 0.01;
    create_output(symnmf(array_of_vectors,H,N,dim_size,300,0.5,eps,num_of_clusters),N,num_of_clusters);
    /*
    if(state == 1){
        create_output(sym(array_of_vectors, N,dim_size),N, N);
    }
    if(state == 2){
        create_output(ddg(array_of_vectors, N, dim_size),N, N);
    }
    if(state == 3){
        create_output(norm(array_of_vectors, N, dim_size),N, N);
    }
     */
    free_matrices(array_of_vectors,N);
    freeVectors(head_vec,N);
    free(head_cord);
    free(curr_vec);
    return 0;

}
