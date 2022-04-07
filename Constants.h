//General
#define DEBUG 0 ///print messages during runtime (doing...done)
#define QUICK 1 ///skips many sorts at expense of more memory usage
#define GIVE_DETAILS 1 ///Prints job info
#define GIVE_ENVS 0 ///prints all unique chemical enviornments

typedef unsigned int uint;

//Input-Output
const extern int BASE;
#define LINESIZE 255 ///max len of any given line in input file

//Size of max element ID string
#define ELEMSIZE 8 ///max number of characters in an element string

//Small element for approximate comparisons
const extern double DRCT_EPS;
const extern double CART_EPS;

//Number of unique enviornments added to the (large) array before 
//we attempt to  re-allocate memory with an additional N_REALLOC
//spaces
const extern unsigned int N_REALLOC;