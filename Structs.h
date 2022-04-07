//Declares Structures and their respective functions
#include "Constants.h"

//A single atomic site
struct site{
	int* elem; ///identifier for this element
	double crdsD[3]; ///direct coordinates a, b, c
	double crdsC[3]; ///cartesian coordinates x, y, z

	struct site* self; ///points to this site's original site
};
double DSqd_s(const struct site* a, const struct site* b);


//A more memory friendly version of a site
//Use when we have massive amounts of sites to hold in memory at once
struct msite{
	short* elem;
	float crds[3];
	struct site *self;
};

//A chemical enviornment comprised of atomic sites
struct env{
	unsigned int nSites;
	struct site** sites;

#if QUICK
	int* sortedElems;
	unsigned int distsLen;
	double* sortedDists;
#endif
};
void MakeElemArr_e(int** arr, 
				   const struct env* e, const unsigned int nElems);
void MakeDstArr_e(double** arr, unsigned int size,
				const struct env* e);
int Cmpr_e(const double* ref, const unsigned int refLen,
		   const struct env* e, const int* elemArr,
		   const unsigned int nElems, const double tol);
int QCmpr_e(const int nElems, const int* refE, const int* tstE,
			const double* ref, const unsigned int refLen,
			const double* tst, const unsigned int tstLen,
			const double tol);
void Write_e(const struct env* ptr);




