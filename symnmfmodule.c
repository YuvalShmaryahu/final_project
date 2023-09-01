#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

static PyObject* fit(PyObject *self, PyObject *args) {
    PyObject *first_matrix, *second_matrix, *current_vector;
    int dim, N, num_of_clusters, state;
    if (!PyArg_ParseTuple(args, "OOii", &first_matrix, &second_matrix, &state,&num_of_clusters)) {
        Py_RETURN_NONE;
    }
    N = PyObject_Length(first_matrix);
    dim = PyObject_Length(PyList_GetItem(first_matrix, 0));
  /* Get the value of the int. */
    int i, j;
    double **X = (double**)calloc(N,sizeof (double*));
    double **H = (double**)calloc(N,sizeof (double*));
    for( i = 0; i < N; i++){
        X[i] = (double *) calloc(dim, sizeof(double));
        current_vector = PyList_GetItem(first_matrix, i);
        for(j = 0; j < dim; j++){
            X[i][j] = PyFloat_AsDouble(PyList_GetItem(current_vector, j));
        }
    }
    for( i = 0; i < N; i++){
        H[i] = (double *) calloc(num_of_clusters, sizeof(double));
        current_vector = PyList_GetItem(second_matrix, i);
        for(j = 0; j < num_of_clusters; j++){
            H[i][j] = PyFloat_AsDouble(PyList_GetItem(current_vector, j));
        }
    }

    if (state == 1){
        H = symnmf(X,H,N,dim,300,0.5,0.0001,num_of_clusters);
    }
    else if (state == 2){
        X = sym(X,N,dim);
    }
    else if (state == 3){
        X = ddg(X,N,dim);
    }
    else {
        X = norm(X,N,dim);
    }
    PyObject* result = PyList_New(N);
    if (state == 1){
        for (i = 0; i < N; i++){
            PyObject* lst = PyList_New(num_of_clusters);
            for(j = 0; j < num_of_clusters; j++){
                PyList_SetItem(lst, j, Py_BuildValue("d", H[i][j]));
            }
            PyList_SetItem(result, i, lst);
        }
    }
    else {
        for (i = 0; i < N; i++){
            PyObject* lst = PyList_New(N);
            for(j = 0; j < N; j++){
                PyList_SetItem(lst, j, Py_BuildValue("d", X[i][j]));
            }
            PyList_SetItem(result, i, lst);
        }
    }
    return result;
}

static PyMethodDef symnmfMethods[]={
        {"fit",
                (PyCFunction) fit,
                   METH_VARARGS,
                     PyDoc_STR("The method expects the following arguments in the following order:\n\n\
                                \t^X=matrixX\n\
                            \tin which each cell in the array is an array of doubles\n\
                            \n\
                            \t^H=matrixH\n\
                            \tin which each cell is an array of doubles\n\
                            \n\
                             t^state = detetmine which function to execute\n\
                            \n\
                            \t^num = num of clusters required")
        },
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef symnmfModule = {
        PyModuleDef_HEAD_INIT,
        "mysymnmfsp",
        NULL,
        -1,
        symnmfMethods
};

PyMODINIT_FUNC
PyInit_mysymnmfsp(void)
{
    PyObject *m;
    m=PyModule_Create(&symnmfModule);
    if(!m)
    {
        Py_RETURN_NONE;
    }
    return m;
}
