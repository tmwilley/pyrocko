

#ifndef PARSTACK_DEFS
#define PARSTACK_DEFS

#include <stdlib.h>

int parstack_config(
        size_t narrays,
        int *offsets,
        size_t *lengths,
        size_t nshifts,
        int *shifts,
        double *weights,
        int method,
        size_t *lengthout,
        int *offsetout);

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
        int nparallel);

#endif
