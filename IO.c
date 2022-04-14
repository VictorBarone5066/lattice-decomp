#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Structs.h"
#include "Decs.h"
#include "Constants.h"

//Sets lattice parameters, element map, site map, cutoff radius
void ReadInputFile(const char* fileName,
                   double ***A,
                   double *rCut,
                   unsigned int *nElems, char ***elems,
                   unsigned int *nSites, struct site ***sites){

	FILE* infile;
	infile = fopen(fileName, "r");
	if (!infile) {
		printf("ReadInputFile(): Unable to locate file :(\n");
		exit(1);
	}


	char line[LINESIZE];
	char* next;

	//Bools for whether we've read in all necessary data
	//Reading loop exits when their sum = 4
	unsigned short rLat = 0; 
	unsigned short rElm = 0; 
	unsigned short rSit = 0;
	unsigned short rRCt = 0;
	for(;(rLat + rElm + rSit + rRCt) != 4;){
		fgets(line, LINESIZE, infile);

		///Read cube domain
		if(strstr(line, "BEGIN LATTICE")){
			*A = malloc(3*sizeof(double*));
			for (unsigned int i = 0; i < 3; ++i) {
				(*A)[i] = malloc(3*sizeof(double));
			}

			for(unsigned int i = 0; i < 3; ++i) {
				fgets(line, LINESIZE, infile);
				(*A)[i][0] = strtod(line, &next);
				(*A)[i][1] = strtod(next, &next);
				(*A)[i][2] = strtod(next, NULL);
			}

			rLat = 1;
		}


		///Read Element Map
		if(strstr(line, "BEGIN ELEMENTMAP")){
			fgets(line, LINESIZE, infile);
			*nElems = strtol(line, &next, BASE);

			*elems = malloc((*nElems)*sizeof(char*));
			for(unsigned int i = 0; i < (*nElems); ++i){
				(*elems)[i] = malloc(ELEMSIZE*sizeof(char));
				fgets(line, LINESIZE, infile);
				for(unsigned int j = 0; j < ELEMSIZE; ++j){
					(*elems)[i][j] = line[j];
				}
			}

			rElm = 1;
		}


		//Read all sites
		if(strstr(line, "BEGIN SITEMAP")){
			fgets(line, LINESIZE, infile);
			*nSites = strtol(line, NULL, BASE);
			*sites = malloc((*nSites)*sizeof(struct site*));
			for(unsigned int i = 0; i < (*nSites); ++i){
				fgets(line, LINESIZE, infile);

				struct site* thisSite = malloc(sizeof(struct site));
				thisSite->elem = malloc(sizeof(int));
				thisSite->self = thisSite;

				thisSite->crdsD[0] = strtod(line, &next);
				thisSite->crdsD[1] = strtod(next, &next);
				thisSite->crdsD[2] = strtod(next, NULL);

				(*sites)[i] = thisSite;
			}

			rSit = 1;
		}

		//Read cutoff radius
		if(strstr(line, "BEGIN CUTOFF")){
			fgets(line, LINESIZE, infile);
			*rCut = strtod(line, NULL);

			rRCt = 1;
		}

	}

	fclose(infile);
	return;
}

//Gives job details at end of run
void PrintDetails(const double **A, const double cutRad,
				  const unsigned long nConfigs,
				  const uint nEnvs, const uint nSites, 
				  const uint nElems,
				  const double runtime){
	printf("\n--------------------------------------------------"
	       "---------------------\n");
	printf("Real-space lattice: (a1x a1y a1z) = "
		   "(%010.5f %010.5f %010.5f)\n", A[0][0], A[0][1], A[0][2]);
	printf("Real-space lattice: (a2x a2y a2z) = "
		   "(%010.5f %010.5f %010.5f)\n", A[1][0], A[1][1], A[1][2]);
	printf("Real-space lattice: (a3x a3y a3z) = "
		   "(%010.5f %010.5f %010.5f)\n", A[2][0], A[2][1], A[2][2]);
	printf("Number of sites: %i\n", nSites);
	printf("Number of elements: %i\n", nElems);
	printf("Cutoff radius: %010.5f\n", cutRad);

	printf("\nTotal runtime: %.1f s\n", runtime);
	printf("Matrix Dimensions = (num configs) x (num unique envs): "
		   "%i x %i\n", nConfigs, nEnvs);
#if WRITE_SPARSE
	printf("Matrix format: sparse (rcvcvcv...)\n");
#endif
#if (!WRITE_SPARSE)
	printf("Matrix format: dense (rvvvvvv...)\n");
#endif
	printf("--------------------------------------------------"
		   "---------------------\n");
}