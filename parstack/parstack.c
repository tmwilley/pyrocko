#include "parstack.h"
#include <omp.h>
#include <stdio.h>

#define CHUNKSIZE 10

int min(int a, int b) {
    return (a < b) ? a : b;
}

int max(int a, int b) {
    return (a > b) ? a : b;
}

double dmax(double a, double b) {
    return (a > b) ? a : b;
}

#define SUCCESS 0
#define NODATA 1

int parstack_config(
        size_t narrays,
        int *offsets,
        size_t *lengths,
        size_t nshifts,
        int *shifts,
        double *weights,
        int method,
        size_t *lengthout,
        int *offsetout) {

    if (narrays < 1) {
        return NODATA;
    }

    int imin, imax, istart, iend;
    size_t iarray, ishift;

    imin = offsets[0] + shifts[0];
    imax = imin + lengths[0];
    for (iarray=0; iarray<narrays; iarray++) {
        for (ishift=0; ishift<nshifts; ishift++) {
            istart = offsets[iarray] + shifts[ishift*narrays + iarray];
            iend = istart + lengths[iarray];
            imin = min(imin, istart);
            imax = max(imax, iend);
        }
    }

    *lengthout = imax - imin;
    *offsetout = imin;

    return SUCCESS;
}

int parstack(
        size_t narrays,
        double **arrays,
        int *offsets,
        size_t *lengths,
        size_t nshifts,
        int *shifts,
        double *weights,
        int method,
        size_t lengthout,
        int offsetout,
        double *result,
        int nparallel) {

    if (narrays < 1) {
        return NODATA;
    }

    int imin, istart;
    size_t iarray, ishift, nsamp, i;
    double weight;
    int chunk;
    double *temp;
    double m;

    imin = offsetout;
    nsamp = lengthout;

    chunk = CHUNKSIZE;

    if (method == 0) {
        #pragma omp parallel private(ishift, iarray, i, istart, weight) num_threads(nparallel)
        {

        #pragma omp for schedule(dynamic,chunk) nowait
        for (ishift=0; ishift<nshifts; ishift++) {
            for (iarray=0; iarray<narrays; iarray++) {
                istart = offsets[iarray] + shifts[ishift*narrays + iarray];
                weight = weights[ishift*narrays + iarray];
                for (i=(size_t)max(0, imin - istart); i<(size_t)max(0, min(nsamp - istart + imin, lengths[iarray])); i++) {
                    result[ishift*nsamp + istart-imin+i] += arrays[iarray][i] * weight;
                }
            }
        }
        }

    } else if (method == 1) {

        #pragma omp parallel private(ishift, iarray, i, istart, weight, temp, m)
        {
        temp = (double*)calloc(nsamp, sizeof(double));
        #pragma omp for schedule(dynamic,chunk) nowait
        for (ishift=0; ishift<nshifts; ishift++) {
            for (i=0; i<nsamp; i++) {
                temp[i] = 0.0;
            }
            for (iarray=0; iarray<narrays; iarray++) {
                istart = offsets[iarray] + shifts[ishift*narrays + iarray];
                weight = weights[ishift*narrays + iarray];
                for (i=(size_t)max(0, imin - istart); i<(size_t)max(0, min(nsamp - istart + imin, lengths[iarray])); i++) {
                    temp[istart-imin+i] += arrays[iarray][i] * weight;
                }
            }
            m = 0.;
            for (i=0; i<nsamp; i++) {
                m += temp[i]*temp[i];
            }
            result[ishift] = m;
        }
        free(temp);
        }
    }
    return SUCCESS;
}

