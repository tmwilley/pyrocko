#include "parstack.h"

int main() {
    size_t narrays = 10;
    double **arrays;
    int *offsets;
    size_t *lengths;
    size_t nshifts;
    int *shifts;
    double *weights;
    double *result;

    size_t iarray, nsamp, ishift, lengthout;
    int offsetout;

    arrays = (double**)calloc(narrays, sizeof(double*));
    offsets = (int*)calloc(narrays, sizeof(int));
    lengths = (size_t*)calloc(narrays, sizeof(size_t));
    for (iarray=0; iarray<narrays; iarray++) {
        nsamp = iarray + 100;
        arrays[iarray] = (double*)calloc(nsamp, sizeof(double));
        lengths[iarray] = nsamp;
        offsets[iarray] = iarray;
    }

    nshifts = 1000000;
    shifts = (int*)calloc(nshifts*narrays, sizeof(int));
    weights = (double*)calloc(nshifts*narrays, sizeof(double));
    for (ishift=0; ishift<nshifts; ishift++) {
        for (iarray=0; iarray<narrays; iarray++) {
            shifts[ishift*narrays +  iarray] = iarray + ishift % narrays;
            weights[ishift*narrays +  iarray] = 1.0;
        }
    }
    parstack_config(narrays, offsets, lengths, nshifts, shifts, weights,
             1, &lengthout, &offsetout);

    result = (double*)calloc(nshifts, sizeof(double));
    parstack(narrays, arrays, offsets, lengths, nshifts, shifts, weights,
             1, lengthout, offsetout, result);


    return 0;
}

