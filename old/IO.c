#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "Decs.h"
#include "Constants.h"
#include "Structs.h"

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


	char* line = malloc(sizeof(char)*LINESIZE);
	char* lineCpy = line;

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
				(*A)[i][0] = strtod(line, &line);
				(*A)[i][1] = strtod(line, &line);
				(*A)[i][2] = strtod(line, NULL);
			}

			rLat = 1;
		}


		///Read Element Map
		if(strstr(line, "BEGIN ELEMENTMAP")){
			fgets(line, LINESIZE, infile);
			*nElems = strtol(line, &line, BASE);

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

				struct site* thisSite; Set_s(&thisSite);
				thisSite->crdsD[0] = strtod(line, &line);
				thisSite->crdsD[1] = strtod(line, &line);
				thisSite->crdsD[2] = strtod(line, NULL);
				
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

	free(lineCpy);
	fclose(infile);
	return;
}