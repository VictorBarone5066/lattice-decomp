//Declares General functions


//File Input/output (IO.c)
void ReadInputFile(const char* fileName,
	               double*** A, double *rCut,
	               unsigned int* nElems, char*** elems,
	               unsigned int* nSites, struct site*** sites);
void PrintDetails(const double** A, const double cutRad,
                  const uint nEnvs, const uint nSites,
                  const uint nElems,
                  const double runtime);

//Related to lattice structure (Lattice.c)
void Inv_3x3(const double** A, double*** AInv);
void SetCarCrds(struct site** s, const double** A);
void SetDirCrds(struct site** s, const double** AInv);
void GiveNImgs(unsigned int** nImgs,
               const double rCut, const double** A);
void MoveToCell(struct site** site);
void SetSiteGeoms(struct env ***envs, const int nEnvs,
                  const unsigned int nSites, 
                  const struct site** sites,
                  const unsigned int* nImgs, 
                  const unsigned int nElems,
                  const double rCut,
                  const double **A);
void SetElems(struct site*** sites, const unsigned int nSites,
              const int* elemArr, const unsigned int elemArrSize);

//Memory management (Memory.c)
void Free_nxm_d(double*** arr,
                const unsigned int n, const unsigned int m);
void Free_nxm_c(char*** arr,
                const unsigned int n, const unsigned int m);