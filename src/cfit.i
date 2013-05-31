%module cfit
/* Include in the wrap code */
%{
#define SWIG_FILE_WITH_INIT
#include "cfit.h"
	%}

%ignore check_input;

/* Numpy magic to output arrays */
%include "numpy.i"
%init %{
	import_array();
	%}

//%apply (int DIM1, double* INPLACE_ARRAY1) {(int gts2, double* frequencies)};
%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int ntp1, int ngt1, int* data)};
%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int ntp2, int ngt2, int* samplesizes)};
%apply (int DIM1, double* IN_ARRAY1) {(int nlocifit, double* fitness)};
%apply (int DIM1, double* IN_ARRAY1) {(int nlociseed, double* seedtimes)};
%apply (int DIM1, double* IN_ARRAY1) {(int nmut, double* mu)};
%apply (int DIM1, int* IN_ARRAY1) {(int ntp, int* tp)};
//%apply (int DIM1, int DIM2, int* OUT_ARRAY2) {(int ntp1, int ngt1, int* data)};

/* Tell SWIG about the declarations */

%include "cfit.h"

