/* -------------------------------------------------------------------
 * FST_DSYMM
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*A*B + BETA*C
 *
 * Exactly one of A, B must be symmetric, in the sense that the UPLO
 * field of the str. code string is specified.
 * NOTE: If C has str. codes, they are ignored.
 *
 * Input:
 * - C:     Input/output matrix. Must have correct size, even if
 *          BETA==0!
 * - A:     Input matrix A
 * - B:     Input matrix B
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

char errMsg[200];

/* Main function FST_DSYMM */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double alpha,beta;
  fst_matrix amat,bmat,cmat;
  bool symm_a;
  char side[2];

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  parseBLASMatrix(prhs[1],"A",&amat,cmat.m,-1);
  parseBLASMatrix(prhs[2],"B",&bmat,amat.n,cmat.n);
  if (UPLO(amat.strcode)!=' ') {
    symm_a=true;
    if (UPLO(bmat.strcode)!=' ')
      mexErrMsgTxt("Exactly one of A, B must have UPLO str. code");
    if (DIAG(amat.strcode)=='U')
      mexErrMsgTxt("DIAG str. code must not be used");
  } else if (UPLO(bmat.strcode)!=' ') {
    symm_a=false;
    if (DIAG(bmat.strcode)=='U')
      mexErrMsgTxt("DIAG str. code must not be used");
  } else
    mexErrMsgTxt("Exactly one of A, B must have UPLO str. code");
  alpha=1.0; beta=0.0;
  if (nrhs>3) {
    alpha=getScalar(prhs[3],"ALPHA");
    if (nrhs>4)
      beta=getScalar(prhs[4],"BETA");
  }

  side[1]=0;
  if (symm_a) {
    side[0]='L';
    BLASFUNC(dsymm) (side,&UPLO(amat.strcode),&cmat.m,&cmat.n,&alpha,
		     amat.buff,&amat.stride,bmat.buff,&bmat.stride,&beta,
		     cmat.buff,&cmat.stride);
  } else {
    side[0]='R';
    BLASFUNC(dsymm) (side,&UPLO(bmat.strcode),&cmat.m,&cmat.n,&alpha,bmat.buff,
		     &bmat.stride,amat.buff,&amat.stride,&beta,cmat.buff,
		     &cmat.stride);
  }
}
