
#define NPY_NO_DEPRECATED_API 7

#include "Python.h"
#include "numpy/arrayobject.h"
#include "parstack.h"

static PyObject *ParstackError;

int good_array(PyObject* o, int typenum) {
    if (!PyArray_Check(o)) {
        PyErr_SetString(ParstackError, "not a NumPy array" );
        return 0;
    }

    if (PyArray_TYPE((PyArrayObject*)o) != typenum) {
        PyErr_SetString(ParstackError, "array of unexpected type");
        return 0;
    }

    if (!PyArray_ISCARRAY((PyArrayObject*)o)) {
        PyErr_SetString(ParstackError, "array is not contiguous or not well behaved");
        return 0;
    }

    return 1;
}

static PyObject* w_parstack(PyObject *dummy, PyObject *args) {

    PyObject *arrays, *offsets, *shifts, *weights, *arr;
    PyObject *result;
    int method, nparallel;
    size_t narrays, nshifts, nweights;
    size_t *clengths;
    size_t lengthout;
    int offsetout;
    int lengthout_arg;
    int *coffsets, *cshifts;
    double *cweights, *cresult;
    double **carrays;
    npy_intp array_dims[1];
    size_t i;
    int err;

    (void)dummy; /* silence warning */

    carrays = NULL;
    clengths = NULL;

    if (!PyArg_ParseTuple(args, "OOOOiiiOi", &arrays, &offsets, &shifts,
                          &weights, &method, &lengthout_arg, &offsetout, &result, &nparallel)) {

        PyErr_SetString(
            ParstackError,
            "usage parstack(arrays, offsets, shifts, weights, method, lengthout, offsetout, result, nparallel)" );

        return NULL;
    }
    if (!good_array(offsets, NPY_INT)) return NULL;
    if (!good_array(shifts, NPY_INT)) return NULL;
    if (!good_array(weights, NPY_DOUBLE)) return NULL;
    if (result != Py_None && !good_array(result, NPY_DOUBLE)) return NULL;

    coffsets = PyArray_DATA(offsets);
    narrays = PyArray_SIZE(offsets);

    cshifts = PyArray_DATA(shifts);
    nshifts = PyArray_SIZE(shifts);

    cweights = PyArray_DATA(weights);
    nweights = PyArray_SIZE(weights);

    nshifts /= narrays;
    nweights /= narrays;

    if (nshifts != nweights) {
        PyErr_SetString(ParstackError, "weights.size != shifts.size" );
        return NULL;
    }

    if (!PyList_Check(arrays)) {
        PyErr_SetString(ParstackError, "arg #1 must be a list of NumPy arrays.");
        return NULL;
    }

    if ((size_t)PyList_Size(arrays) != narrays) {
        PyErr_SetString(ParstackError, "len(offsets) != len(arrays)");
        return NULL;
    }

    carrays = (double**)calloc(narrays, sizeof(double*));
    if (carrays == NULL) {
        PyErr_SetString(ParstackError, "alloc failed");
        return NULL;
    }

    clengths = (size_t*)calloc(narrays, sizeof(size_t));
    if (clengths == NULL) {
        PyErr_SetString(ParstackError, "alloc failed");
        free(carrays);
        return NULL;
    }

    for (i=0; i<narrays; i++) {
        arr = PyList_GetItem(arrays, i);
        if (!good_array(arr, NPY_DOUBLE)) {
            free(carrays);
            free(clengths);
            return NULL;
        }
        carrays[i] = PyArray_DATA(arr);
        clengths[i] = PyArray_SIZE(arr);
    }
    if (lengthout_arg < 0) {
        err = parstack_config(narrays, coffsets, clengths, nshifts, cshifts,
                              cweights, method, &lengthout, &offsetout);

        if (err != 0) {
            PyErr_SetString(ParstackError, "parstack_config() failed");
            free(carrays);
            free(clengths);
            return NULL;
        }
    } else {
        lengthout = (size_t)lengthout_arg;
    }

    if (method == 0) {
        array_dims[0] = nshifts * lengthout;
    } else {
        array_dims[0] = nshifts;
    }

    if (result != Py_None) {
        if (PyArray_SIZE(result) != array_dims[0]) {
            free(carrays);
            free(clengths);
            return NULL;
        }
        Py_INCREF(result); 
    } else {
        result = PyArray_SimpleNew(1, array_dims, NPY_FLOAT64);
        cresult = PyArray_DATA(result);

        for (i=0; i<(size_t)array_dims[0]; i++) {
            cresult[i] = 0.0;
        }

        if (result == NULL) {
            free(carrays);
            free(clengths);
            return NULL;
        }
    }
    cresult = PyArray_DATA(result);

    err = parstack(narrays, carrays, coffsets, clengths, nshifts, cshifts,
                   cweights, method, lengthout, offsetout, cresult, nparallel);

    if (err != 0) {
        PyErr_SetString(ParstackError, "parstack() failed");
        free(carrays);
        free(clengths);
        Py_DECREF(result);
        return NULL;
    }

    free(carrays);
    free(clengths);
    return Py_BuildValue("Ni", result, offsetout);
}

static PyMethodDef ParstackMethods[] = {
    {"parstack",  w_parstack, METH_VARARGS,
        "Parallel weight-and-delay stacking" },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC initparstack_ext(void) {
    PyObject *m;

    m = Py_InitModule("parstack_ext", ParstackMethods);
    if (m == NULL) return;
    import_array();

    ParstackError = PyErr_NewException("parstack_ext.error", NULL, NULL);
    Py_INCREF(ParstackError);  /* required, because other code could remove `error` 
                               from the module, what would create a dangling
                               pointer. */
    PyModule_AddObject(m, "ParstackError", ParstackError);
}

