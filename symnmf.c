#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
                A[i][j] = exp((-euclidean_distance(X[i], X[j], N)) / 2);
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
