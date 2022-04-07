//Declares Structures and their respective functions
#include "Constants.h"

//A single atomic site
struct site{
	int* elem; ///identifier for this element
	double crdsD[3]; ///direct coordinates a, b, c
	double crdsC[3]; ///cartesian coordinates x, y, z
	
	//unsigned int nNbrs;
	//struct site **nbrs; ///array of this site's neighbors

	struct site *self; ///points to this site's original site
};
void Set_s(struct site** ptr);
void DeepCopy_s(struct site** sNew, const struct site *sOld);
void Free_s(struct site** ptr);

//A chemical enviornment comprised of atomic sites
struct env{
	unsigned int nSites;
	struct site** sites;
};
void Set_e(struct env** ptr,
           const unsigned int nNbrs);
void SetRefEnv(struct env** ptr, const struct env* ref,
			   const unsigned int nElems);
int Cmpr_e(const struct env* a, const struct env* b,
		   const double tol);
void Free_e(struct env** ptr);

