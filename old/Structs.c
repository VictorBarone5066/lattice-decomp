//Implements Structures and their respective functions
#include <stdlib.h>
#include <math.h>

#include "Structs.h"
#include "Constants.h"
#include "SpecLinAlg.h"

//A single atomic site
void Set_s(struct site** ptr){
	*ptr = malloc(sizeof(struct site));
	//(*ptr)->nNbrs = 0;
	//(*ptr)->nbrs = malloc(0);
	(*ptr)->self = *ptr;
	
	return;
}

void DeepCopy_s(struct site** sNew, const struct site* sOld){
	*sNew = malloc(sizeof(struct site));
	(*sNew)->elem = sOld->elem;
	for(unsigned int i = 0; i < 3; ++i){
		(*sNew)->crdsC[i] = sOld->crdsC[i];
		(*sNew)->crdsD[i] = sOld->crdsD[i];
	}

	//(*sNew)->nNbrs = sOld->nNbrs;
	//(*sNew)->nbrs = malloc(((*sNew)->nNbrs)*sizeof(struct site));
	//for(unsigned int i = 0; i < (*sNew)->nNbrs; ++i){
	//	(*sNew)->nbrs[i] = sOld->nbrs[i];	
	//}

	(*sNew)->self = sOld->self;

	return;
}

void Free_s(struct site** ptr){
	//free((*ptr)->nbrs);
	free((*ptr));

	return;
}


//A chemical enviornment comprised of atomic sites
void Set_e(struct env** ptr,
           const unsigned int nSites){
	*ptr = malloc(sizeof(struct env));
	(*ptr)->nSites = nSites;
	(*ptr)->sites = malloc(nSites*sizeof(struct site*));
}

void DeepCopy_e(struct env** eNew, const struct env* eOld){
	Set_e(&(*eNew), eOld->nSites);
	for(unsigned int i = 0; i < (*eNew)->nSites; ++i){
		struct site* tmp;
		DeepCopy_s(&tmp, eOld->sites[i]);
		(*eNew)->sites[i] = tmp;
	}
}


//Returns the 4D COM as cent{centX, centY, centZ, centE}
//cent[4] should be allocated to 0.0 
void COM_e(double** cent, const struct env* self){
	for(unsigned int i = 0; i < self->nSites; ++i){
		for(unsigned int j = 0; j < 3; ++j){
			(*cent)[j] += self->sites[i]->crdsC[j];
		}
		(*cent)[3] += (double)(*self->sites[i]->elem);
	}
	for(unsigned int i = 0; i < 4; ++i){
		(*cent)[i] /= (double)(self->nSites);
	}
}

//Compares two enviornments by TODO: FINISH DESCRIPTION WHEN I FIGURE OUT HOW I WANT THIS TO WORK
int Cmpr_e(const struct env* a, const struct env* b, 
			  const double tol){
	///Setup data matrices
	///A is 4 x nSites, and B is nSites x 4 so that we do C = A^T B
	///Do the shifts to the respective COMS as well
	double** AT = malloc(4*sizeof(double*));
	double** B = malloc((b->nSites)*sizeof(double*));
	double* comA = calloc(4, sizeof(double));
	double* comB = calloc(4, sizeof(double));
	COM_e(&comA, a); COM_e(&comB, b);

	for(unsigned int i = 0; i < 3; ++i){
		AT[i] = malloc((a->nSites)*sizeof(double));
		for(unsigned int j = 0; j < (a->nSites); ++j){
			AT[i][j] = a->sites[j]->crdsC[i] - comA[i];
		}
	}
	AT[3] = malloc((a->nSites)*sizeof(double));
	for(unsigned int i = 0; i < (a->nSites); ++i){
		AT[3][i] = (double)(*a->sites[i]->elem) - comA[3];
	}

	for(unsigned int i = 0; i < (b->nSites); ++i){
		B[i] = malloc(4*sizeof(double));
		for(unsigned int j = 0; j < 3; ++j){
			B[i][j] = b->sites[i]->crdsC[j] - comB[j];
		}
		B[i][3] = (double)(*b->sites[i]->elem) - comB[3];
	}

	///Need to do SVD on the correlation of A and B
	double** C; 
	IProd_2(&C, AT, 4, a->nSites, B, b->nSites, 4);

	///Compute best rotation
	/// !!! NOTE !!! It would be fairly easy to optimize this: The 
	/// cross correlation matrix C (and therefore all SVD matrices
	/// related to it) will be 4 x 4 - just do some #defines and such
	/// Also, it would be useful to make SVD() return U^T instead of 
	/// U.  Also Also, we can save some inverse sqrts in SVD() since
	/// we never actually use the S matrix.  
	double** U; double** S; double** V;
	SVD(&U, &S, &V, C, 4, 4, 5);
	double** UT;
	Transpose(&UT, U, 4, 4);
	double** R;
	IProd_2(&R, V, 4, 4, UT, 4, 4);

	double** BRot;
	//IProd_2(&BRot, B, b->nSites, 4, R, 4, 4);


	///The rotation matrix should be a 4x4 identity if the two envs
	///are equal.  
	///Additionally, the rotation matrix is guarenteed to be 
	///normalized, so if |the sum of the main diagonals| =/= 4, the
	///two envs can not be equal
	int ret = 0;
	double dev = R[0][0] + R[1][1] + R[2][2] + R[3][3];
	if(fabs(dev - 4.0) < tol){
		ret = 1;
	}

	///Clean up
	free(comA);
	free(comB);
	Free_nxm_d(&AT, 4, a->nSites);
	Free_nxm_d(&B, b->nSites, 4);
	Free_nxm_d(&C, 4, 4);
	Free_nxm_d(&U, 4, 4);
	Free_nxm_d(&UT, 4, 4);
	Free_nxm_d(&R, 4, 4);

	return ret;
}

void SetRefEnv(struct env** ptr, const struct env* ref,
			   const unsigned int nElems){
	DeepCopy_e(&(*ptr), ref);
	for(unsigned int i = 0; i < ref->nSites; ++i){
		*(*ptr)->sites[i]->elem = (i)%nElems;
	}
}

void Free_e(struct env** ptr){
	for(unsigned int i = 0; i < (*ptr)->nSites; ++i){
		Free_s(&((*ptr)->sites[i]));
	}
	free((*ptr)->sites);
	free((*ptr));
}