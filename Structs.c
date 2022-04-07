//Implements Structures and their respective functions
#include <stdlib.h>
#include <math.h>

#include "Structs.h"
#include "Constants.h"

int dblcmpr(const void* a, const void* b);

//A single atomic site
void Malloc_s(struct site** ptr){
	*ptr = malloc(sizeof(struct site));

	return;
}

//Squared distance between two sites
//!!! Treats the element index as the 4th dimension !!!
//It may later become necessary that each element pair (a, b) maps to
//a unique value c.  If (a, b) -> c and (b, a) -> c is acceptable, 
//then the function c = (a+b) + (a*b) works well for small (a, b).
double DSqd_s(const struct site* a, const struct site* b){
	return (
		    (b->crdsC[0] - a->crdsC[0])*(b->crdsC[0] - a->crdsC[0]) +
		    (b->crdsC[1] - a->crdsC[1])*(b->crdsC[1] - a->crdsC[1]) +
		    (b->crdsC[2] - a->crdsC[2])*(b->crdsC[2] - a->crdsC[2]) +
			((double)*b->elem - (double)*a->elem)*
								((double)*b->elem - (double)*a->elem)
		   );
}

//Makes an array where the ith component provides the count of how 
//many elements of type i are in an enviornment
void MakeElemArr_e(int** arr, 
				   const struct env* e, const unsigned int nElems){
	*arr = calloc(nElems, sizeof(int));
	for(unsigned int i = 0; i < e->nSites; ++i){
		(*arr)[*(e->sites[i]->elem)] += 1;
	}
}

//Makes a sorted distance array from an enviornment,
//fills the *arr pointer with it
void MakeDstArr_e(double** arr, const unsigned int size,
				  const struct env* e){
	*arr = malloc(size*sizeof(double));
	for(unsigned int i = 0, loc = 0; i < e->nSites - 1; ++i){
		for(unsigned int j = i + 1; j < e->nSites; ++j, ++loc){
			(*arr)[loc] = DSqd_s(e->sites[i], e->sites[j]);
		}
	}
	qsort(*arr, size, sizeof(double), dblcmpr);
}


//Compares two enviornments
//A set of distances between all points in a body uniquely defines 
//a body, up to some translation and rotation (the proof of this is
//left to the reader).
//So, two {vectors} are chemically equivalent if their sorted 
//distances all match up.
//Takes an already sorted array of distances (the reference) and 
//compares it to an enviornment by finding the distances in the 
//enviornment, sorting it, and comparing the two arrays
int dblcmpr(const void* a, const void* b){
	if(*(const double*)a < *(const double*)b) return -1;
	if(*(const double*)a > *(const double*)b) return +1;
	return 0;
}
// !!! TODO: I have a feeling that this will break for arrays of
//           size < 3 or 4 ... double check that all bases are 
//           covered !!!
int Cmpr_e(const double* ref, const unsigned int refLen,
		   const struct env* e, const int* elemArr,
		   const unsigned int nElems, const double tol){
	unsigned int tstSize = (e->nSites)*(e->nSites - 1)/2;
	///Don't even bother if the arrays are differently sized
	if(refLen != tstSize){
		return 0;
	}
	///Don't even bother if the arrays have different numbers of
	///elements
	int* tstE;
	MakeElemArr_e(&tstE, e, nElems);
	for(unsigned int i = 0; i < nElems; ++i){
		if(elemArr[i] != tstE[i]){
			free(tstE);
			return 0;
		}
	}

	///Make array for this enviornment, sort it
	double *tstD;
	MakeDstArr_e(&tstD, tstSize, e);

	///Compare element-by-element
	for(unsigned int i = 0; i < tstSize; ++i){
		if(fabs(ref[i] - tstD[i]) > tol){
			free(tstD);
			free(tstE);
			return 0;
		}
	}

	free(tstD);
	free(tstE);
	return 1;
}

//Quickly compare envs (don't sort every time...)
//Assumes ref and tstLen are already sorted. 
int QCmpr_e(const int nElems, const int* refE, const int* tstE,
			const double* ref, const unsigned int refLen,
			const double* tst, const unsigned int tstLen,
			const double tol){
	if(refLen != tstLen) {
		return 0;
	}

	for(unsigned int i = 0; i < nElems; ++i) {
		if (refE[i] != tstE[i]) {
			return 0;
		}
	}

	for(unsigned int i = 0; i < refLen; ++i) {
		if (fabs(ref[i] - tst[i]) > tol) {
			return 0;
		}
	}

	return 1;
}

void Write_e(const struct env* ptr){
	for(unsigned int i = 0; i < ptr->nSites; ++i){
		printf("Elem: %i | Crds(a, b, c): (%f, %f, %f)\n",
			*ptr->sites[i]->elem,
			ptr->sites[i]->crdsD[0],
			ptr->sites[i]->crdsD[1],
			ptr->sites[i]->crdsD[2]);
	}
}




