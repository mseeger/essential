/* -------------------------------------------------------------------
 * FST_DSYRK
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*op(A)*op(A)' + BETA*C, op(X) == X or X', dep. on TRS
 *
 * Here, C must be symmetric (UPLO str. code required). A str. codes
 * are ignored.
 *
 * Input:
 * - C:     Input/output matrix. Must have correct size. Must be
 *          symmetric (UPLO code req.)
 * - A:     Input matrix A
 * - TRS:   Flag, s.a.
 * - ALPHA: S.a. Def.: 1
 * - BETA:  S.a. Def.: 0
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include "mex.h"
#include "mex_helper.h"
#include "blas_headers.h"

char errMsg[200];

/* Main function FST_DGEMM */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int k;
  double alpha,beta;
  fst_matrix amat,cmat;
  char trans[2];

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  if (UPLO(cmat.strcode)==' ')
    mexErrMsgTxt("C must be symmetric (UPLO code req.)");
  if (DIAG(cmat.strcode)=='U')
    mexErrMsgTxt("DIAG str. code must not be used");
  trans[1]=0;
  if (getScalInt(prhs[2],"TRS")!=0) {
    trans[0]='T';
    parseBLASMatrix(prhs[1],"A",&amat,-1,cmat.n);
    k=amat.m;
  } else {
    trans[0]='N';
    parseBLASMatrix(prhs[1],"A",&amat,cmat.n,-1);
    k=amat.n;
  }
  alpha=1.0; beta=0.0;
  if (nrhs>3) {
    alpha=getScalar(prhs[3],"ALPHA");
    if (nrhs>4)
      beta=getScalar(prhs[4],"BETA");
  }

  BLASFUNC(dsyrk) (&UPLO(cmat.strcode),trans,&cmat.n,&k,&alpha,amat.buff,
		   &amat.stride,&beta,cmat.buff,&cmat.stride);
}
