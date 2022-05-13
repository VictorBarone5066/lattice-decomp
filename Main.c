#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Structs.h"
#include "Decs.h";
#include "Constants.h"



char INFILE_NAME[LINESIZE] = "input.ldc";
double RAD_OVERWRITE = -1.0;

int main(int argc, char *argv[]){
#if GIVE_DETAILS
	clock_t clockStart = clock();
#endif

	//CLIs-----------------------------------------------------------
	for(unsigned int i = 0; i < argc; ++i){
		//Infile name
		if((argv[i][0] == '-') && 
		   (argv[i][1] == 'i' || argv[i][1] == 'I'))
			for(unsigned int j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					break;
				}
				INFILE_NAME[j] = argv[i + 1][j];
			}
		//Overwritten cutoff radius
		if((argv[i][0] == '-') &&
		   (argv[i][1] == 'r' || argv[i][1] == 'R')){
			RAD_OVERWRITE = atof(argv[i + 1]);
		}

			
	}
	//---------------------------------------------------------------


	//Initial Input--------------------------------------------------
	double** A; double rCut;
	uint nElems; char** elems;
	uint nSites; struct site** sites;
#if DEBUG
	printf("reading infile ... ");
#endif
	ReadInputFile(INFILE_NAME, &A, &rCut, 
		          &nElems, &elems, 
		          &nSites, &sites);
	if(RAD_OVERWRITE >= +0.0){
		rCut = RAD_OVERWRITE;
	}
#if DEBUG
	printf("done\n");
#endif
	//---------------------------------------------------------------


	//Assign site geometries-----------------------------------------
	uint* nImgs; GiveNImgs(&nImgs, rCut, A);
	for(uint i = 0; i < nSites; ++i){MoveToCell(&sites[i]);}
	
	struct env** envs;
#if DEBUG
	printf("setting sites ... ");
#endif
	SetSiteGeoms(&envs, nSites, nSites, sites, 
				 nImgs, nElems, rCut, A);
	if(N_REALLOC < nSites){
		printf("main(): Increase N_REALLOC and re-compile\n");
		exit(1);
	}
#if DEBUG
	printf("done\n");
#endif
	//---------------------------------------------------------------


	//Setup done: Now loop through all configs-----------------------
	unsigned long nConfigs;
#if DEBUG
	printf("finding unique envs ... \n");
#endif
	FILE* infile;
	infile = fopen(INFILE_NAME, "r");
	if (!infile) {
		printf("main(): Unable to locate file :(\n");
		exit(1);
	}
	char line[LINESIZE];
	char* next;

	///Stuff for keeping track of unique enviornments
	unsigned long uniqueCurSize = 0;
	unsigned long uniqueMaxSize = N_REALLOC;
	struct env** uniqueEnvs = malloc(N_REALLOC*sizeof(struct env*));

	///Stuff for writing out each cell's enviornment decomposition
	///cellDecomp[a] = b means that this cell has b instances of
	///enviornment number a
	uint* cellDecomp = malloc(N_REALLOC*sizeof(uint));


	for(;;){
	fgets(line, LINESIZE, infile);
	if(strstr(line, "BEGIN CONFIG")){
	fgets(line, LINESIZE, infile);
	nConfigs = strtoul(line, &next, BASE);
	int* elemArr = malloc(nSites*sizeof(int));

	for(unsigned long i = 0; i < nConfigs; ++i){
		///Fill element array
		fgets(line, LINESIZE, infile);
		next = line;
		for(uint j = 0; j < nSites; ++j){
			elemArr[j] = strtol(next, &next, BASE);
		}

		///Prep decomposition vector
		memset(cellDecomp, 0, uniqueMaxSize*sizeof(uint));

		///Set sites with new elements...
		SetElems(&sites, nSites, elemArr, nSites);

		///...and check for any unique enviornments
		for(uint j = 0; j < nSites; ++j){
			////comparison arrays to avoid re-making env[j]'s
			////arrays every time we look at uniqueEnv[k]
			////(for speed)
			double* refArr; int* elemArr_;
			int refArrSize = (envs[j]->nSites)*
								(envs[j]->nSites - 1)/2;
			MakeDstArr_e(&refArr, refArrSize, envs[j]);
			MakeElemArr_e(&elemArr_, envs[j], nElems);

			///Check all known unique envs to see if the current env
			///exists there or not
			unsigned short toAdd = 1;
			for(unsigned long k = 0; k < uniqueCurSize; ++k){
#if QUICK
				if(QCmpr_e(nElems, elemArr_,
					       uniqueEnvs[k]->sortedElems, refArr,
					       refArrSize, uniqueEnvs[k]->sortedDists,
					       uniqueEnvs[k]->distsLen, CART_EPS)){
					toAdd = 0;
					cellDecomp[k]++;
					break;
			}
#endif
#if (!QUICK)
				if(Cmpr_e(refArr, refArrSize, uniqueEnvs[k], 
						  elemArr_, nElems, CART_EPS)){
					toAdd = 0;
					cellDecomp[k]++;
					break;
				}
#endif
			}

			//If no identical envs were found in the current array...
			if(toAdd){
				cellDecomp[uniqueCurSize]++;

				///Make deep copy of a new enviornment to add to 
				///the unique env array
				struct env* newEnv = malloc(sizeof(struct env));
				newEnv->nSites = envs[j]->nSites;
				newEnv->sites = malloc(envs[j]->nSites*
										sizeof(struct site));
				for(uint k = 0; k < newEnv->nSites; ++k){
					///Deepcopy each individual site (yikes)
					struct site* tSite = malloc(sizeof(struct site));
					tSite->elem = malloc(sizeof(int));
					*tSite->elem = *envs[j]->sites[k]->elem;
					for(uint m = 0; m < 3; ++m){
						tSite->crdsD[m]= envs[j]->sites[k]->crdsD[m];
						tSite->crdsC[m]= envs[j]->sites[k]->crdsC[m];
					}
					tSite->self = envs[j]->sites[k]->self;

					newEnv->sites[k] = tSite;
				}
#if QUICK
				int* elmsarr;
				MakeElemArr_e(&elmsarr, newEnv, nElems);
				newEnv->sortedElems = elmsarr;
				newEnv->distsLen = (newEnv->nSites)*(newEnv->nSites - 1)/2;
				double *darr;
				MakeDstArr_e(&darr, newEnv->distsLen, newEnv);
				newEnv->sortedDists = darr;
#endif
				uniqueEnvs[uniqueCurSize] = newEnv;
				uniqueCurSize++;
			}

			free(refArr);
			free(elemArr_);
		}
				
		///Print this config's env vector
		printf("%i", i);
		for(uint j = 0; j < uniqueCurSize; ++j){
#if WRITE_SPARSE
			if(cellDecomp[j] != 0){
				printf(" %i %i", j, cellDecomp[j]);
			}
#endif
#if (!WRITE_SPARSE)
			printf(" %i", cellDecomp[j]);
#endif
		}
		printf("\n");

		///Make more room for enviornments if necessary
		if(uniqueMaxSize - uniqueCurSize < nSites){
			uniqueMaxSize += N_REALLOC;
			uniqueEnvs = realloc(uniqueEnvs, 
							 	 uniqueMaxSize*sizeof(struct env*));
			cellDecomp = realloc(cellDecomp, 
							     uniqueMaxSize*sizeof(uint));
		}
	}

	//we've finished reading all of the configs now
	free(elemArr);
	goto BreakRead;
	}
	}

	BreakRead:
	fclose(infile);
#if DEBUG
	printf("done\n");
#endif
	//---------------------------------------------------------------

#if GIVE_ENVS
	printf("\n--------------------------------------------------"
		"---------------------\n");
	for(uint i = 0; i < uniqueCurSize; ++i){
		Write_e(uniqueEnvs[i]);
		printf("\n");
	}
	printf("--------------------------------------------------"
		"---------------------\n");
#endif

	//Clean up-------------------------------------------------------
#if DEBUG
	printf("cleaning up ... ");
#endif
	for(uint i = 0; i < uniqueCurSize; ++i){
		for(uint j = 0; j < uniqueEnvs[i]->nSites; ++j){
			free(uniqueEnvs[i]->sites[j]->elem);
			free(uniqueEnvs[i]->sites[j]);
		}
#if QUICK
		free(uniqueEnvs[i]->sortedElems);
		free(uniqueEnvs[i]->sortedDists);
#endif
		free(uniqueEnvs[i]->sites);
		free(uniqueEnvs[i]);
	}
	free(uniqueEnvs);
	for(uint i = 0; i < nSites; ++i){
		for(uint j = 0; j < envs[i]->nSites; ++j){
			///no, I havn't forgotten to free elem - that's later
			free(envs[i]->sites[j]);
		}
		free(envs[i]->sites);
		free(envs[i]);

		free(sites[i]->elem);
		free(sites[i]);
	}
	free(cellDecomp);
	free(envs);
	free(sites);
	Free_nxm_c(&elems, nElems, ELEMSIZE);
	free(nImgs);
#if GIVE_DETAILS
	clock_t clockStop = clock();
	double runtime = (double)(clockStop - clockStart)/CLOCKS_PER_SEC;
	PrintDetails(A, rCut, nConfigs, uniqueCurSize, nSites, nElems, runtime);
#endif
	Free_nxm_d(&A, 3, 3);
#if DEBUG
	printf("done\n");
#endif

	return 0;
}