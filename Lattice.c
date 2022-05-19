//Functions for lattice related work
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Structs.h"
#include "Decs.h"

#ifdef __unix__
#include "extern/kdtree-master/kdtree.h"
#endif
#ifdef _WIN32
#include "extern//kdtree-master//kdtree.h"
#endif


//Inverse of the 3x3 lattice matrix A
//A = { {ax, ay, az}, {bx, by, bz}, {cx, cy, cz} }
//where a, b, c are the first, second, and third lattice vectors
//The inverse is set s.t. the transform to direct is done in the 
//same way that the transform to cartesian is done.  
double Det_2x2(const double tl, const double tr,
	           const double bl, const double br) {
	return tl*br - tr*bl;
}
void Inv_3x3(const double** A, double*** AInv){
	*AInv = malloc(3 * sizeof(double*));
	for (unsigned int i = 0; i < 3; ++i) {
		(*AInv)[i] = malloc(3 * sizeof(double));
	}

	double det = A[0][0]*Det_2x2(A[1][1], A[1][2], A[2][1], A[2][2])
               - A[0][1]*Det_2x2(A[1][0], A[1][2], A[2][0], A[2][2])
               + A[0][2]*Det_2x2(A[1][0], A[1][1], A[2][0], A[2][1]);
	///no, I didn't mess this up - im applying the negative signs 
	///and transpose all in one step
	(*AInv)[0][0] = Det_2x2(A[1][1], A[1][2], A[2][1], A[2][2])/det;
	(*AInv)[0][1] = -Det_2x2(A[0][1], A[0][2], A[2][1], A[2][2])/det;
	(*AInv)[0][2] = Det_2x2(A[0][1], A[0][2], A[1][1], A[1][2])/det;
	(*AInv)[1][0] = -Det_2x2(A[1][0], A[1][2], A[2][0], A[2][2])/det;
	(*AInv)[1][1] = Det_2x2(A[0][0], A[0][2], A[2][0], A[2][2])/det;
	(*AInv)[1][2] = -Det_2x2(A[0][0], A[0][2], A[1][0], A[1][2])/det;
	(*AInv)[2][0] = Det_2x2(A[1][0], A[1][1], A[2][0], A[2][1])/det;
	(*AInv)[2][1] = -Det_2x2(A[0][0], A[0][1], A[2][0], A[2][1])/det;
	(*AInv)[2][2] = Det_2x2(A[0][0], A[0][1], A[1][0], A[1][1])/det;

	return;
}

//Set Cartesian coordinates given direct coordinates and the lattice
void SetCarCrds(struct site** s, const double** A){
	for(unsigned int i = 0; i < 3; ++i){
		(*s)->crdsC[i] = 0.0;
		for(unsigned int j = 0; j < 3; ++j){
			(*s)->crdsC[i] += (*s)->crdsD[j]*A[j][i]; 
		}
	}

	return;
}

//Set Direct coordinates given cartesian coordinates and the inv latt
void SetDirCrds(struct site** s, const double** AInv){
	for(unsigned int i = 0; i < 3; ++i) {
		(*s)->crdsD[i] = 0.0;
		for(unsigned int j = 0; j < 3; ++j) {
			(*s)->crdsD[i] += (*s)->crdsC[j]*AInv[j][i];
		}
	}

	return;
}

//Returns the number of images necessary for periodic BCs.  
//We need enough periodic images s.t. rCut at any edge does not pass
//a periodic wall.  
void GiveNImgs(unsigned int** nImgs, 
               const double rCut, const double **A){
	*nImgs = malloc(3*sizeof(int));
	for(unsigned int i = 0; i < 3; ++i){
		double latLen = sqrt(A[i][0]*A[i][0] + 
                             A[i][1]*A[i][1] + 
                             A[i][2]*A[i][2]);
		for(unsigned int j = 1;;++j){
			if(rCut < (double)j*latLen){
				(*nImgs)[i] = j;
				break;
			}
		}
	}
}

//Changes direct coordinates of a site s.t. 0 <= site coord i < 1
void MoveToCell(struct site **site){
	for(unsigned short i = 0; i < 3; ++i){
		while((*site)->crdsD[i] < 0.0){
			(*site)->crdsD[i] += 1.0;
		}
		while((*site)->crdsD[i] >= 1.0){
			(*site)->crdsD[i] -= 1.0;
		}
	}

	return;
}


void SetSiteGeoms(struct env*** envs, const int nEnvs,
				  const unsigned int nSites, 
                  const struct site** sites,
				  const unsigned int *nImgs, 
				  const unsigned int nElems,
                  const double rCut,
                  const double **A){
	///Get array of all neighbors 
	unsigned int addSize = (2*nImgs[0] + 1)*
						   (2*nImgs[1] + 1)*
						   (2*nImgs[2] + 1)*
						   nSites;
	struct site** additionals = malloc(addSize*sizeof(struct site*));
	unsigned int count = 0;
	for(int a = -((int)nImgs[0]); a <= (int)nImgs[0]; ++a){
	for(int b = -((int)nImgs[1]); b <= (int)nImgs[1]; ++b){
	for(int c = -((int)nImgs[2]); c <= (int)nImgs[2]; ++c){
		for(unsigned int i = 0; i < nSites; ++i){
			struct site* thisSite = malloc(sizeof(struct site)); 
			for (unsigned int j = 0; j < 3; ++j) {
				thisSite->crdsD[j] = sites[i]->crdsD[j];
				thisSite->crdsC[j] = sites[i]->crdsC[j];
			}
			thisSite->elem = sites[i]->elem;
			thisSite->self = sites[i]->self;


			thisSite->crdsD[0] += (double)a;
			thisSite->crdsD[1] += (double)b;
			thisSite->crdsD[2] += (double)c;
			SetCarCrds(&thisSite, A);

			///-1: original, -2: periodic image
			thisSite->elem = (a == 0 && b == 0 && c == 0)? -1:-2;

			additionals[count] = thisSite;
			count++;
		}
	}
	}
	}

	//Connect sites to periodic images with a KD tree
	void* kdNode = kd_create(3);
	for(unsigned int i = 0; i < addSize; ++i){
		kd_insert3(kdNode,
				   additionals[i]->crdsC[0],
				   additionals[i]->crdsC[1],
				   additionals[i]->crdsC[2],
				   i);
	}
	
	///Set the list of enviornments
	*envs = malloc(nEnvs*sizeof(struct env*));
	double x = 0.0; double y = 0.0; double z = 0.0;
	unsigned int envCount = 0;
	for(unsigned int i = 0; i < addSize; ++i){
		////Assign nbrs only to sites that lie within the orig cell
		if (additionals[i]->elem != -1) {
			continue;
		}
		////temp variables - kd node CHANGES THE COORDINATES
		////somewhere in their functions...
		x = additionals[i]->crdsC[0];
		y = additionals[i]->crdsC[1];
		z = additionals[i]->crdsC[2];

		struct kdres* res = kd_nearest_range3(kdNode, x, y, z, rCut);
		int resSize = kd_res_size(res);

		//Set_e(&thisEnv, resSize);
		struct env* thisEnv = malloc(sizeof(struct env));
		thisEnv->nSites = resSize;
		thisEnv->sites = malloc(resSize*sizeof(struct site*));

		////Assign neighbors to this enviornment
		for(int j = 0; j < resSize; ++j){
			unsigned int thisInd = (unsigned int*)kd_res_item3(res, &x, &y, &z);
			
			struct site* nbr = malloc(sizeof(struct site));
			nbr->elem = additionals[thisInd]->elem;
			for(unsigned int k = 0; k < 3; ++k){
				nbr->crdsD[k] = additionals[thisInd]->crdsD[k];
				nbr->crdsC[k] = additionals[thisInd]->crdsC[k];
			}
			nbr->self = additionals[thisInd]->self;
			nbr->elem = nbr->self->elem;

			thisEnv->sites[j] = nbr;
			kd_res_next(res);
		}

		thisEnv->nSites = resSize;
		(*envs)[envCount] = thisEnv;

		kd_res_free(res);
		envCount++;
	}


	kd_free(kdNode);
	for(unsigned int i = 0; i < addSize; ++i){
		free(additionals[i]);
	}
	free(additionals);
	
	return;
}

//Takes an array of ints (elemArr) of size elemArrSize to set the
//elements of the enviornment array (through editing the sites, not
//the env. array - the env array's elements point to sites in the 
//unit cell so this works fine)
void SetElems(struct site*** sites, const unsigned int nSites,
			  const int* elemArr, const unsigned int elemArrSize){
	//The element at site i is equal to the integer at the ith 
	//component in the element array
	for(unsigned int i = 0; i < elemArrSize; ++i){
		*(*sites)[i]->elem = elemArr[i];
	}
}