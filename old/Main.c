#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Decs.h";
#include "Constants.h"
#include "Structs.h"
#include "SpecLinAlg.h"



char INFILE_NAME[LINESIZE] = "input";


int main(int argc, char *argv[]){
	
	//Initial Input--------------------------------------------------
	double** A; double rCut;
	unsigned int nElems; char** elems;
	unsigned int nSites; struct site** sites;
	ReadInputFile(INFILE_NAME, &A, &rCut, 
		          &nElems, &elems, 
		          &nSites, &sites);
	//---------------------------------------------------------------


	//Assign site geometries-----------------------------------------
	double** AInv; Inv_3x3(A, &AInv);
	unsigned int* nImgs; GiveNImgs(&nImgs, rCut, A);
	for(unsigned int i = 0; i < nSites; ++i){MoveToCell(&sites[i]);}
	
	struct env** envs; struct env* ref;
	SetSiteGeoms(&envs, nSites, nSites, sites, 
				 nImgs, nElems, rCut, A);
	SetRefEnv(&ref, envs[0], nElems);
	//---------------------------------------------------------------


	//Setup done: Now loop through all configs-----------------------
	FILE* infile;
	infile = fopen(INFILE_NAME, "r");
	if (!infile) {
		printf("main(): Unable to locate file :(\n");
		exit(1);
	}
	char* line = malloc(sizeof(char)*LINESIZE);
	char* lineCpy = line;

	for(;;){
		fgets(line, LINESIZE, infile);
		if(strstr(line, "BEGIN CONFIG")){
			fgets(line, LINESIZE, infile);
			unsigned long nConfigs = strtoul(line, &line, BASE);
			int* elemArr = malloc(nSites*sizeof(int));

			for(unsigned long i = 0; i < nConfigs; ++i){
				///Fill element array
				fgets(line, LINESIZE, infile);
				for(unsigned int j = 0; j < nSites; ++j){
					elemArr[j] = strtol(line, &line, BASE);
				}

				///Set sites with new elements
				SetElems(&sites, nSites, elemArr, nSites);
				printf(":)");
			}

			//we've finished reading all of the configs now
			free(elemArr);
			goto BreakRead;
		}
	}

	BreakRead:
	free(lineCpy);
	fclose(infile);
	//---------------------------------------------------------------

	int dumb;
	dumb = Cmpr_e(ref, envs[0], 0.01);
	struct site* tmp;
	tmp = envs[0]->sites[0];
	envs[0]->sites[0] = envs[0]->sites[1];
	envs[0]->sites[1] = tmp;
	dumb = Cmpr_e(ref, envs[0], 0.01);




	//
	int y = 0;
	sites[2]->elem = 3;
	//Clean up
	printf("end");
	Free_e(&ref);
	Free_nxm_c(&elems, nElems, ELEMSIZE);
	Free_nxm_e(&envs, nSites, 1); 
	Free_nxm_s(&sites, nSites, 1);
	Free_nxm_d(&A, 3, 3);
	Free_nxm_d(&AInv, 3, 3);
	free(nImgs);

	return 0;
}