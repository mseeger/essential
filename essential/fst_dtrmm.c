/* -------------------------------------------------------------------
 * FST_DTRMM
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * B = ALPHA*op(A)*B (LEFT==1), op(X) == X or X', dep. on TRS
 * B = ALPHA*B*op(A) (LEFT==0)
 *
 * Here, A mu best triangular (UPLO str. code required, DIAG code
 * optional) and square. If B has str. codes, they are ignored.
 *
 * Input:
 * - B:     Input/output matrix
 * - A:     Input matrix A (triangular; UPLO code req.)
 * - LEFT:  Flag, s.a.
 * - TRS:   Flag, s.a.
 * - ALPHA: S.a. Def.: 1
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

/* Main function FST_DTRMM */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  double alpha;
  fst_matrix amat,bmat;
  char side[2],trans[2];

  /* Read arguments */
  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"B",&bmat,-1,-1);
  side[1]=0;
  if (getScalInt(prhs[2],"LEFT")!=0) {
    side[0]='L';
    parseBLASMatrix(prhs[1],"A",&amat,bmat.m,bmat.m);
  } else {
    side[0]='R';
    parseBLASMatrix(prhs[1],"A",&amat,bmat.n,bmat.n);
  }
  if (UPLO(amat.strcode)==' ')
    mexErrMsgTxt("A must be triangular (UPLO code req.)");
  trans[1]=0;
  trans[0]=(getScalInt(prhs[3],"TRS")!=0)?'T':'N';
  alpha=1.0;
  if (nrhs>4)
    alpha=getScalar(prhs[4],"ALPHA");

  BLASFUNC(dtrmm) (side,&UPLO(amat.strcode),trans,&DIAG(amat.strcode),&bmat.m,
		   &bmat.n,&alpha,amat.buff,&amat.stride,bmat.buff,
		   &bmat.stride);
}
