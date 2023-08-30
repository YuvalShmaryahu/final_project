#include <Python.h>
#include "symnmf.h"

static PyObject *wrap_function(PyObject *self, PyObject *args, double** (*func)(double**, int, int)) {
    PyObject *lst;
    PyObject *curr_vec;
    PyObject *answer;
    double **vectors_mat;
    int vec_num, vec_len, i, j, p;

    if (!PyArg_ParseTuple(args, "iiO", &vec_len, &vec_num, &lst)) {
        return NULL;
    }

    vectors_mat = create_vectors_mat(vec_num, vec_len);
    if (vectors_mat == NULL) {
        return NULL;
    }

    for (j = 0; j < vec_num; ++j) {
        curr_vec = PyList_GetItem(lst, j);
        for (p = 0; p < vec_len; ++p) {
            vectors_mat[j][p] = PyFloat_AsDouble(PyList_GetItem(curr_vec, p));
        }
    }

    double **result_mat = func(vectors_mat, vec_num, vec_len);

    answer = PyList_New(vec_num);
    for (i = 0; i < vec_num; ++i) {
        curr_vec = PyList_New(vec_len);
        for (j = 0; j < vec_len; ++j) {
            PyList_SetItem(curr_vec, j, Py_BuildValue("d", result_mat[i][j]));
        }
        PyList_SetItem(answer, i, curr_vec);
    }

    free_mat(result_mat, vec_num);
    free_mat(vectors_mat, vec_num);
    return answer;
}

static PyObject *py_symnmf(PyObject *self, PyObject *args) {
    return wrap_function(self, args, finalmat);
}

static PyObject *py_sym(PyObject *self, PyObject *args) {
    return wrap_function(self, args, symmat);
}

static PyObject *py_ddg(PyObject *self, PyObject *args) {
    return wrap_function(self, args, ddgmat);
}

static PyObject *py_norm(PyObject *self, PyObject *args) {
    return wrap_function(self, args, normmat);
}

static PyMethodDef symnmfMethods[] = {
        {"symnmf", (PyCFunction) py_symnmf, METH_VARARGS, PyDoc_STR("function expects: vector's dimension, number of vectors, initialized H, W")},
        {"sym", (PyCFunction) py_sym, METH_VARARGS, PyDoc_STR("function expects: vector's dimension, number of vectors, vectors matrix")},
        {"ddg", (PyCFunction) py_ddg, METH_VARARGS, PyDoc_STR("function expects: vector's dimension, number of vectors, vectors matrix")},
        {"norm", (PyCFunction) py_norm, METH_VARARGS, PyDoc_STR("function expects: vector's dimension, number of vectors, vectors matrix")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef Symnmfmodule = {
        PyModuleDef_HEAD_INIT,
        "symnmfmodule",
        NULL,
        -1,
        symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *m;
    m = PyModule_Create(&Symnmfmodule);
    if (!m) {
        Py_RETURN_NONE;
    }
    return m;
}

