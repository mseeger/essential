/* -------------------------------------------------------------------
 * INVCHOL
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * Overwrites Cholesky factor of matrix A by inverse of A. The
 * factor is given in FACT, which is upper or lower triangular
 * (UPLO field required). The inverse overwrites FACT, but just the
 * triangle occupied by the factor is accessed.
 * Uses the LAPACK routine DPOTRI.
 *
 * Input:
 * - FACT: Cholesky factor of A, overwritten by inverse of A
 * Return:
 * - INFO: Return code of DPOTRI (0: success)
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include "mex.h"
#include "mex_helper.h"

char errMsg[200];

/* LAPACK DPOTRI declaration */
int BLASFUNC(dpotri) (const char* uplo,const int* n,double* a,const int* lda,
		      int* info);

/* Main function FST_INVCHOL */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,n;
  fst_matrix fact;
  char uplo[2];

  /* Read arguments */
  if (nrhs<1)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"FACT",&fact,-1,-1);
  if ((n=fact.n)!=fact.m || UPLO(fact.strcode)==' ')
    mexErrMsgTxt("FACT must be triangular");

  /* Call DPOTRI */
  BLASFUNC(dpotri) (&UPLO(fact.strcode),&n,fact.buff,&fact.stride,&i);
  /* Return INFO */
  if (nlhs>0) {
#ifdef MATLAB_VER65
    plhs[0]=mxCreateDoubleScalar((double) i);
#else
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=(double) i;
#endif
  }
}
