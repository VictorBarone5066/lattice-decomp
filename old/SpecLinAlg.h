//Declares Specalized linear algebra functions
#include "Constants.h"

//Basic Stuff
double IProd_1(const int len,
			   const double* u, const double* v);
void IProd_2(double*** res,
			 const double** A, const int aRow, const int aCol,
			 const double** B, const int bRow, const int bCol);
void Transpose(double*** res, const double** M,
			   const int nRow, const int nCol);

//Specialized stuff
void JacobiRight(double*** M, const int size,
				 const double c, const double s,
				 const int p, const int q);
void JacobiRight(double*** M, const int size,
				 const double c, const double s,
				 const int p, const int q);
void JacobiSandwich(double*** M, const int size,
					const double c, const double s,
					const int p, const int q);
void SymSchurDecomp(double* c, double* s,
					const double** M,
					const int p, const int q);
void CyclicJacobi(double*** M, double*** V, 
				  const int size, const int nSweeps);
void SVD(double*** U, double*** S, double*** V,
		 const double** A, const int m, const int n,
		 const int nSweeps);
