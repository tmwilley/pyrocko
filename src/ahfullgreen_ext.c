#define NPY_NO_DEPRECATED_API 7

#include "Python.h"
#include "numpy/arrayobject.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef npy_float64 float64_t;

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

typedef enum {
    SUCCESS = 0,
    SINGULARITY,
    BAD_ARRAY,
} ahfullgreen_error_t;

const char* ahfullgreen_error_names[] = {
    "SUCCESS",
    "SINGULARITY",
    "BAD_ARRAY",
};

static PyObject *Error;


int good_array(PyObject* o, int typenum, ssize_t size_want) {
    if (!PyArray_Check(o)) {
        PyErr_SetString(Error, "not a NumPy array" );
        return 0;
    }

    if (PyArray_TYPE((PyArrayObject*)o) != typenum) {
        PyErr_SetString(Error, "array of unexpected type");
        return 0;
    }

    if (!PyArray_ISCARRAY((PyArrayObject*)o)) {
        PyErr_SetString(Error, "array is not contiguous or not well behaved");
        return 0;
    }

    if (size_want != -1 && size_want != PyArray_SIZE((PyArrayObject*)o)) {
        PyErr_SetString(Error, "array is of wrong size");
        return 0;
    }

    return 1;
}

static ahfullgreen_error_t numpy_or_none_to_c_double(
        PyObject* o, ssize_t size_want, double **arr, size_t *size) {

    if (o == Py_None) {
        if (size_want > 0) {
            PyErr_SetString(Error, "array is of wrong size");
            return BAD_ARRAY;
        }
        *arr = NULL;
        *size = 0;
    } else {
        if (!good_array(o, NPY_DOUBLE, -1)) return BAD_ARRAY;
        *arr = PyArray_DATA((PyArrayObject*)o);
        *size = PyArray_SIZE((PyArrayObject*)o);
    }

    return SUCCESS;
}

static ahfullgreen_error_t numpy_or_none_to_c_complex(
        PyObject* o, ssize_t size_want, double complex **arr, size_t *size) {

    if (o == Py_None) {
        if (size_want > 0) {
            PyErr_SetString(Error, "array is of wrong size");
            return BAD_ARRAY;
        }
        *arr = NULL;
        *size = 0;
    } else {
        if (!good_array(o, NPY_COMPLEX128, -1)) return BAD_ARRAY;
        *arr = PyArray_DATA((PyArrayObject*)o);
        *size = PyArray_SIZE((PyArrayObject*)o);
    }

    return SUCCESS;
}

static ahfullgreen_error_t add_seismogram(
        double vp,
        double vs,
        double density,
        double qp,
        double qs,
        double *x,
        double *f,
        double *m6,
        int out_quantity,  // 1: velocity, 2: acceleration
        double out_delta,
        double out_offset,
        size_t out_size,
        double complex *out_x,
        double complex *out_y,
        double complex *out_z,
        int want_far,
        int want_intermediate,
        int want_near
        ) {

    double r, r2, r4, density4pi;
    double gamma[3];
    double complex *out[3];
    double complex *b1, *b2, *b3;
    double m[3][3] = {{m6[0], m6[3], m6[4]},
         {m6[3], m6[1], m6[5]},
         {m6[4], m6[5], m6[2]}};

    int n, p, q, i;

    double a1, a2, a3, a4, a5, a6, a7, a8;
    double vp2, vp3, vs2, vs3, w;
    double complex iw, dfactor;


    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    if (r == 0.0) {
        return SINGULARITY;
    }

    r2 = r*r;
    r4 = r2*r2;
    vs2 = vs*vs;
    vs3 = vs2*vs;
    vp2 = vp*vp;
    vp3 = vp2*vp;

    density4pi = density * M_PI * 4.;

    for (n=0; n<3; n++) gamma[n] = x[n]/r;

    out[0] = out_x;
    out[1] = out_y;
    out[2] = out_z;

    b1 = (double complex*)calloc(out_size, sizeof(double complex));
    b2 = (double complex*)calloc(out_size, sizeof(double complex));
    b3 = (double complex*)calloc(out_size, sizeof(double complex));

    for (i=0; i<out_size; i++) {
        w = out_offset + out_delta * i;
        iw = I * w;
        dfactor = 1.0/iw;
        if (out_quantity == 1) {
            dfactor = 1.0;
        } else if (out_quantity == 2) {
            dfactor = iw;
        }
        if (i != 0) {
            b2[i] = dfactor * cexp(-iw * r/vp) * exp(-w * r / (2.0*vp*qp));
            b3[i] = dfactor * cexp(-iw * r/vs) * exp(-w * r / (2.0*vs*qs));
            b1[i] = (r/vp + 1.0/iw) * b2[i]/iw - (r/vs + 1.0/iw) * b3[i]/iw;
        } else {
            b2[i] = 0.0;
            b3[i] = 0.0;
            b1[i] = 0.0;
        }
    }

    for (n=0; n<3; n++) {
        if (out[n] == NULL) continue;

        for (p=0; p<3; p++) {
            for (q=0; q<3; q++) {

                if (want_near) {
                    a1 = (
                        15. * gamma[n] * gamma[p] * gamma[q] -
                        3. * (gamma[n] * (p==q) +
                            gamma[p] * (n==q) +
                            gamma[q] * (n==p))) /
                        (density4pi * r4);
                } else {
                    a1 = 0.;
                }

                if (want_intermediate) {
                    a2 = (
                        6. * gamma[n] * gamma[p] * gamma[q] -
                        gamma[n] * (p==q) -
                        gamma[p] * (n==q) -
                        gamma[q] * (n==p)) /
                        (density4pi * vp2 * r2);

                    a3 = - (
                        6. * gamma[n] * gamma[p] * gamma[q] -
                        gamma[n] * (p==q) -
                        gamma[p] * (n==q) -
                        2. * gamma[q] * (n==p)) /
                        (density4pi * vs2 * r2);
                } else {
                    a2 = a3 = 0.;
                }

                if (want_far) {
                    a4 = (gamma[n] * gamma[p] * gamma[q]) /
                        (density4pi * vp3 * r);


                    a5 = - (gamma[q] * (gamma[n] * gamma[p] - (n==p))) /
                        (density4pi * vs3 * r);
                } else {
                    a4 = a5 = 0.;
                }

                for (i=0; i<out_size; i++) {
                    iw = I * (out_offset + out_delta * i);
                    out[n][i] += (a1*b1[i] + a2*b2[i] + a3*b3[i] +
                                iw*a4*b2[i] + iw*a5*b3[i]) * m[p][q];
                }
            }
        }

        for (p=0; p<3; p++) {

            a6 = (3. * gamma[n] * gamma[p] - (n==p)) /
                (density4pi * r2 * r);

            a7 = (gamma[n] * gamma[p]) /
                (density4pi * vp2 * r);

            a8 = - (gamma[n] * gamma[p] - (n==p)) /
                (density4pi * vs2 * r);

            for (i=0; i<out_size; i++) {
                out[n][i] += (a6*b1[i] + a7*b2[i] + a8*b3[i]) * f[p];
            }
        }
    }

    free(b3);
    free(b2);
    free(b1);

    return SUCCESS;
}

static PyObject* w_add_seismogram(PyObject *dummy, PyObject *args) {
    double vp, vs, density, qp, qs;
    double *x;
    double *f;
    double *m6;
    int out_quantity;  // 0: displacement, 1: velocity, 2: acceleration
    double out_delta;
    double out_offset;

    size_t out_size, out_x_size, out_y_size, out_z_size;

    PyObject *x_arr;
    PyObject *f_arr;
    PyObject *m6_arr;

    PyObject *out_x_arr;
    PyObject *out_y_arr;
    PyObject *out_z_arr;

    double complex *out_x;
    double complex *out_y;
    double complex *out_z;

    int want_far;
    int want_intermediate;
    int want_near;

    ahfullgreen_error_t err;
    size_t dummy_size;

    (void)dummy; /* silence warning */

    if (!PyArg_ParseTuple(args, "dddddOOOiddOOOiii",
            &vp, &vs, &density, &qp, &qs, &x_arr, &f_arr, &m6_arr,
            &out_quantity, &out_delta, &out_offset,
            &out_x_arr, &out_y_arr, &out_z_arr, &want_far, &want_intermediate,
            &want_near)) {

        PyErr_SetString(Error,
            "usage: add_seismogram(vp, vs, density, qp, qs, x, f, m6, "
            "out_quantity, out_delta, out_offset, "
            "out_x, out_y, out_z, want_far, want_intermediate, want_near)");

        return NULL;
    }

    if (SUCCESS != numpy_or_none_to_c_complex(out_x_arr, -1, &out_x, &out_x_size)) return NULL;
    if (SUCCESS != numpy_or_none_to_c_complex(out_y_arr, -1, &out_y, &out_y_size)) return NULL;
    if (SUCCESS != numpy_or_none_to_c_complex(out_z_arr, -1, &out_z, &out_z_size)) return NULL;

    out_size = max(max(out_x_size, out_y_size), out_z_size);

    if ((!(out_x_size == 0 || out_x_size == out_size)) ||
        (!(out_y_size == 0 || out_y_size == out_size)) ||
        (!(out_z_size == 0 || out_z_size == out_size))) {

        PyErr_SetString(Error, "differing output array sizes");
        return NULL;
    }

    if (SUCCESS != numpy_or_none_to_c_double(x_arr, 3, &x, &dummy_size)) return NULL;
    if (SUCCESS != numpy_or_none_to_c_double(f_arr, 3, &f, &dummy_size)) return NULL;
    if (SUCCESS != numpy_or_none_to_c_double(m6_arr, 6, &m6, &dummy_size)) return NULL;

    err = add_seismogram(
        vp, vs, density, qp, qs, x, f, m6,
        out_quantity, out_delta, out_offset, out_size,
        out_x, out_y, out_z, want_far, want_intermediate, want_near);

    if (err != SUCCESS) {
        PyErr_SetString(Error, ahfullgreen_error_names[err]);
        return NULL;
    }

    return Py_BuildValue("");
}

static PyMethodDef StoreExtMethods[] = {
    {"add_seismogram", w_add_seismogram, METH_VARARGS,
        "Add seismogram to array." },

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initahfullgreen_ext(void)
{
    PyObject *m;

    m = Py_InitModule("ahfullgreen_ext", StoreExtMethods);
    if (m == NULL) return;
    import_array();

    Error = PyErr_NewException("pyrocko.ahfullgreen_ext.Error", NULL, NULL);
    Py_INCREF(Error);  /* required, because other code could remove `error`
                               from the module, what would create a dangling
                               pointer. */
    PyModule_AddObject(m, "Error", Error);
}

