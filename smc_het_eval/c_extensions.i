%module c_extensions

// Numpy Related Includes:
%{
#define SWIG_FILE_WITH_INIT
#include "c_extensions.h"
%}

// numpy arrays
%include "numpy.i"
%init %{
import_array(); // This is essential. We will get a crash in Python without it.
%}

// interface
%apply (double * INPLACE_ARRAY2,int DIM1, int DIM2 ) {(double * x, int dimx1, int dimx2)};
%apply (double * IN_ARRAY1, int DIM1) {(double * mask, int lenmask)};
%include "c_extensions.h"
