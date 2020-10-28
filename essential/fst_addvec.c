/* -------------------------------------------------------------------
 * FST_ADDVEC
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*X*1' + C  (COL==1)
 * C = ALPHA*1*X' + C  (COL==0)
 *
 * X vector, 1 vector of all ones. Structure codes ignored.
 *
 * Input:
 * - C:     Input/output matrix. Must have correct size, even if
 *          BETA==0!
 * - X:     Input vector X
 * - COL:   S.a.
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

/* Main function FST_ */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,ione=1;
  double alpha;
  fst_matrix cmat;
  const double* xvec;
  bool col;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  col=(getScalInt(prhs[2],"COL")!=0);
  i=getVecLen(prhs[1],"X");
  if ((col && i!=cmat.m) || (!col && i!=cmat.n))
    mexErrMsgTxt("X has wrong size");
  xvec=mxGetPr(prhs[1]);
  alpha=1.0;
  if (nrhs>3)
    alpha=getScalar(prhs[3],"ALPHA");

  if (col) {
    for (i=0; i<cmat.n; i++)
      BLASFUNC(daxpy) (&cmat.m,&alpha,xvec,&ione,cmat.buff+i*cmat.stride,
		       &ione);
  } else {
    for (i=0; i<cmat.m; i++)
      BLASFUNC(daxpy) (&cmat.n,&alpha,xvec,&ione,cmat.buff+i,
		       &cmat.stride);
  }
}
