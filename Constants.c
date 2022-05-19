#include "Constants.h"

typedef unsigned int uint;

//File Input-Output
const extern int BASE = 10;


//Small element for approximate comparisons
const extern double DRCT_EPS = 1.0E-2; ///fractional 
const extern double CART_EPS = 5.0E-1; ///angstroms

//Number of unique enviornments added to the (large) array before 
//we attempt to  re-allocate memory with an additional N_REALLOC
//spaces
//Should be AT LEAST as large as the number of sites per cell
const extern unsigned int N_REALLOC = 256;