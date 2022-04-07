//For memory management.  Mostly freeing stuff
#include <stdlib.h>

#include "Decs.h"
#include "Structs.h"

//Free double array of size n x m
//m is not necessary, but it makes it look neater...
void Free_nxm_d(double*** arr, 
	            const unsigned int n, const unsigned int m){
	for(unsigned int i = 0; i < n; ++i){
		free((*arr)[i]);
	}
	free((*arr));
}

//Free char array of size n x m
//m is not necessary, but it makes it look neater...
void Free_nxm_c(char*** arr,
                const unsigned int n, const unsigned int m) {
	for(unsigned int i = 0; i < n; ++i){
		free((*arr)[i]);
	}
	free((*arr));
}

//Free sites array of size n x m
//m is not necessary, but it makes it look neater...
void Free_nxm_s(struct site*** arr,
                const unsigned int n, const unsigned int m){
	for(unsigned int i = 0; i < n; ++i) {
		Free_s(&(*arr)[i]);
	}
	free((*arr));
}

void Free_nxm_e(struct env ***arr,
				const unsigned int n, const unsigned int m){
	for(unsigned int i = 0; i < n; ++i){
		Free_e(&(*arr)[i]);
	}
	free((*arr));
}