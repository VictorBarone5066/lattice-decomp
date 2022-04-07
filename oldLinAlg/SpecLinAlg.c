//Implements Specalized linear algebra functions
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Decs.h"
#include "Constants.h"
#include "SpecLinAlg.h"

//Inner product of two real 1-d vectors
double IProd_1(const int len,
			   const double *u, const double *v){
	double res = 0.0;
	for(unsigned int i = 0; i < len; ++i){
		res += u[i]*v[i]; 
	}
	return res;
}

//Matrix product of two real 2-d (square) matrices 
//stores the result in 'res': res = AB
//Note: This may very well be a large source of running time and is
//likley due to cache misses instead of actual calculations...
//...However, the sizes of A and B are pretty small in this program
//(about 12).  A more effective product R_ij = sum_k (A_ik*B_kj) 
//could be to first transpose B, do the multiplication
void IProd_2(double ***res, 
			 const double **A, const int aRow, const int aCol, 
			 const double **B, const int bRow, const int bCol){
	*res = malloc(aRow*sizeof(double*));
	for(unsigned int i = 0; i < aRow; ++i){
		(*res)[i] = malloc(bCol*sizeof(double));
	}

	for(unsigned int i = 0; i < aRow; ++i){
		for(unsigned int j = 0; j < bCol; ++j){
			(*res)[i][j] = 0.0;
			for(unsigned int k = 0; k < aCol; ++k){
				(*res)[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return;
}

//Transpose of matrix M stored in res
//nRow, nCol are the rows and columns of M, NOT RES
void Transpose(double ***res, const double **M,
			   const int nRow, const int nCol){
	*res = malloc(nCol*sizeof(double*));
	for (unsigned int i = 0; i < nCol; ++i){
		(*res)[i] = malloc(nRow*sizeof(double));
	}
	
	double tmp = 0.0;
	for(unsigned int i = 0; i < nRow; ++i){
		for(unsigned int j = 0; j < nCol; ++j){
			(*res)[j][i] = M[i][j];
		}
	}

	return;
}


//Functions for diagonalization of a SYMMETRIC matrix M of size n by
//the cyclic jacobi method w/ schur decompositions
//Reference: Matrix Computations, 3rd ed. by Golub & Van Loan

//Transforms M -> JM where J is the jacobi rotation matrix 
//with cosine part c and sine part s at indices p < q
void JacobiLeft(double*** M, const int size,
				const double c, const double s,
				const int p, const int q){
	///make a temp copy of M that won't be changed 
	double** tmp;
	tmp = malloc(size * sizeof(double*));
	for (unsigned int i = 0; i < size; ++i) {
		tmp[i] = malloc(size * sizeof(double));
		memcpy(tmp[i], (*M)[i], size * sizeof(double));
	}

	///Only the rows of M are affected by J
	for(unsigned int j = 0; j < size; ++j){
		(*M)[p][j] = c*tmp[p][j] - s*tmp[q][j];
		(*M)[q][j] = c*tmp[q][j] + s*tmp[p][j];
	}
	
	///clean up
	Free_nxm_d(&tmp, size, size);
	return;
}

//Transforms M -> MJ where J is the jacobi rotation matrix 
//with cosine part c and sine part s at indices p < q
void JacobiRight(double ***M, const int size,
			     const double c, const double s,
				 const int p, const int q){
	///make a temp copy of M that won't be changed 
	double** tmp;
	tmp = malloc(size*sizeof(double*));
	for (unsigned int i = 0; i < size; ++i) {
		tmp[i] = malloc(size*sizeof(double));
		memcpy(tmp[i], (*M)[i], size*sizeof(double));
	}
	

	///Only the columns of M are affected by J
	for(unsigned int i = 0; i < size; ++i){
		(*M)[i][p] = c*tmp[i][p] + s*tmp[i][q];
		(*M)[i][q] = c*tmp[i][q] - s*tmp[i][p];
	}

	///clean up
	Free_nxm_d(&tmp, size, size);
	return;
}

//Transforms M -> J^t M J where J are the jacobi rotation matrices 
//with cosine part c and sine part s at indices p < q
void JacobiSandwich(double ***M, const int size,
					const double c, const double s,
					const int p, const int q){
	///make a temp copy of M that won't be changed 
	double** tmp;
	tmp = malloc(size*sizeof(double*));
	for(unsigned int i = 0; i < size; ++i){
		tmp[i] = malloc(size*sizeof(double));
		memcpy(tmp[i], (*M)[i], size*sizeof(double));
	}


	///Loop over all M[p/q][i] = M[i][p/q] before setting specific 
	///M[p][q] elements to avoid annoying if-else chains
	for(unsigned int i = 0; i < size; ++i){
		(*M)[p][i] = c*tmp[p][i] - s*tmp[q][i];
		(*M)[i][p] = (*M)[p][i];

		(*M)[q][i] = s*tmp[p][i] + c*tmp[q][i];
		(*M)[i][q] = (*M)[q][i];
	}
	///Now, specific M[p][q] elements
	double cc = c*c;
	double ss = s*s;
	double cs = c*s;
	(*M)[p][p] = cc*tmp[p][p] - 2.0*cs*tmp[p][q] + ss*tmp[q][q];
	(*M)[q][q] = ss*tmp[p][p] + 2.0*cs*tmp[p][q] + cc*tmp[q][q];
	(*M)[p][q] = (cc-ss)*tmp[p][q] + cs*(tmp[p][p]-tmp[q][q]);
	(*M)[q][p] = (*M)[p][q];

	///clean up
	Free_nxm_d(&tmp, size, size);
	return;
}

//Compute the cosine-sine pair (c, s) from Schur Decomp given matrix
//indices 1 <= p < q <= n
void SymSchurDecomp(double *c, double *s, 
					const double** M,
					const int p, const int q){
	if(M[p][q] == 0.0){
		*c = 1.0;
		*s = 0.0;
		return;
	}

	double tau = (M[q][q] - M[p][p])/(2.0*M[p][q]);
	double mult = (tau >= 0.0)? +1.0: -1.0;
	double t = mult/(mult*tau + sqrt(1.0 + tau*tau));

	*c = 1.0/sqrt(1.0 + t*t);
	*s = t*(*c);
	return;
}

//Reduces M to "almost diagonal" through the cyclic Jacobi algo.  
//Will preform nSweeps passes (5 is usually safe)
//Stores accumulated jacobi rotations in matrices V (right)
void CyclicJacobi(double*** M, double*** V, 
				  const int size, const int nSweeps){
	///Initialize an ortogonal vector V to the identity
	// !!! NOTE: Potential Disaster !!!
	//C does not guarentee that the memory calloced to doubles will
	//represent zero - this is system-dependent.  In most cases 
	//(i.e. the ones I care about), it ends up being zero, though.  
	*V = malloc(size*sizeof(double*));
	for(unsigned int i = 0; i < size; ++i){
		(*V)[i] = calloc(size, sizeof(double));
		(*V)[i][i] = 1.0;
	}
	
	double c; double s;
	for(unsigned int n = 0; n < nSweeps; ++n){
		for(int p = 0; p < size - 1; ++p){
			for(int q = p + 1; q < size; ++q){
				SymSchurDecomp(&c, &s, (*M), p, q);
				JacobiLeft(&(*V), size, c, s, p, q);
				JacobiSandwich(&(*M), size, c, s, p, q);
			}
		}
	}

	return;
}

//Helper function to SVD() that sorts two ragged matrices in place 
//according to the eigenvalues
//Eigenvalues stored in E, eigenvectors stored in V
void DoubleSort(double*** E, double*** V, const int n){
	for (int i = 0; i < n - 1; ++i) {
		int a = i;
		for (int j = i + 1; j < n; ++j) {
			if ((*E)[j][j] > (*E)[a][a]) {
				a = j;
			}
		}
		if (i != a) {
			///Swap eigenvalues a and i
			double tmpEig = (*E)[a][a];
			(*E)[a][a] = (*E)[i][i];
			(*E)[i][i] = tmpEig;

			///Swap eigenvectors a and i (note: the use of ragged 
			///arrays here makes this very easy) :)
			double* tmpRow = (*V)[a];
			(*V)[a] = (*V)[i];
			(*V)[i] = tmpRow;
		}
	}
}

//Returns the SVD(A_nxm) -> U, S, V using the Jacobi method
//with nSweeps number of passes
//Note U's, V's ROWS are their eigenvectors - V is NOT transposed 
//If A is size mxn, then U = mxm, S = mxn, V = nxn
void SVD(double ***U, double ***S, double ***V, 
		 const double **A, const int m, const int n,
		 const int nSweeps){
	double** AT;
	Transpose(&AT, A, m, n);

	///*******************************
	///CASE num rows(m) > num cols (n)
	///*******************************
	if(m > n){
		//First, Compute A A^T  
		double** AAT;
		IProd_2(&AAT, A, m, n, AT, n, m);
		//and diagonalize it, accumulating eigenvectors as rows in V
		CyclicJacobi(&AAT, &(*U), m, nSweeps);

		//Sort to follow SVD formalism
		DoubleSort(&AAT, &(*U), m);

		//S is diagonal with the ith eigenvalue of A A^T
		//on the ith diagonal spot
		*S = malloc(m*sizeof(double*));
		for(unsigned int i = 0; i < m; ++i){
			(*S)[i] = calloc(n, sizeof(double));
			(*S)[i][i] = sqrt(fabs(AAT[i][i]));
		}

		//The ith row of V is equal to 
		// 1/sig_i * A^T dot (the ith eigenvector of A A^T) which are just
		//the rows of U
		*V = malloc(n*sizeof(double*));
		for(unsigned int i = 0; i < n; ++i){
			///Dot product bc for some reason IProd_2() dosn't work here
			(*V)[i] = calloc(n, sizeof(double));
			if((*S)[i][i] < SINGVAL_TOL){
				continue;
			}
			for(unsigned int j = 0; j < n; ++j){
				for(unsigned int k = 0; k < m; ++k){
					(*V)[i][j] += AT[j][k] * (*U)[i][k];
				}
				///scale
				(*V)[i][j] *= 1.0/((*S)[i][i]);
			}
		}

		Free_nxm_d(&AT, n, m); 
		Free_nxm_d(&AAT, m, m);
		return;
	}





	///*******************************
	///CASE otherwise (i.e. m <= n)
	///*******************************
	//First, Compute A^T A 
	double** ATA;
	IProd_2(&ATA, AT, n, m, A, m, n);
	//and diagonalize it, accumulating eigenvectors as rows in V
	CyclicJacobi(&ATA, &(*V), n, nSweeps);

	//Sort to follow SVD formalism
	DoubleSort(&ATA, &(*V), n);

	//S is diagonal with the ith eigenvalue of A^T A
	//on the ith diagonal spot
	*S = malloc(m*sizeof(double*));
	for(unsigned int i = 0; i < m; ++i){
		(*S)[i] = calloc(n, sizeof(double));
		(*S)[i][i] = sqrt(fabs(ATA[i][i]));
	}

	//The ith row of U is equal to 
	// 1/sig_i * A dot (the ith eigenvector of A^T A) which are just
	//the rows of V
	*U = malloc(m*sizeof(double*));
	for(unsigned int i = 0; i < m; ++i){
		///Dot product bc for some reason IProd_2() dosn't work here
		(*U)[i] = calloc(m, sizeof(double));
		if((*S)[i][i] < SINGVAL_TOL){
			continue;
		}
		for(unsigned int j = 0; j < m; ++j){
			for(unsigned int k = 0; k < n; ++k){
				(*U)[i][j] += A[j][k] * (*V)[i][k];	
			}
			///scale
			(*U)[i][j] *= 1.0/((*S)[i][i]);
		}
	}

	Free_nxm_d(&AT, n, m);
	Free_nxm_d(&ATA, n, n);
	return;
}
