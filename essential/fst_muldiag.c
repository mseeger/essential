/* -------------------------------------------------------------------
 * FST_MULDIAG
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * A = DIAG(D)*A (LEFT==1)
 * A = A*DIAG(D) (LEFT==0)
 *
 * Input:
 * - A:    Input/output matrix
 * - D:    Input vector (diag. matrix)
 * - LEFT: Flag (s.a.)
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

/* Main function FST_MULDIAG */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,ione=1;
  const double* dvec;
  fst_matrix amat;
  bool left;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"A",&amat,-1,-1);
  if (amat.m*amat.n==0)
    mexErrMsgTxt("Wrong argument A");
  left=(getScalInt(prhs[2],"LEFT")!=0);
  i=getVecLen(prhs[1],"D");
  if ((left && i!=amat.m) || (!left && i!=amat.n))
    mexErrMsgTxt("Arguments have wrong sizes");
  dvec=mxGetPr(prhs[1]);

  if (left) {
    for (i=0; i<amat.m; i++)
      BLASFUNC(dscal) (&amat.n,dvec+i,amat.buff+i,&amat.stride);
  } else {
    for (i=0; i<amat.n; i++)
      BLASFUNC(dscal) (&amat.m,dvec+i,amat.buff+(i*amat.stride),&ione);
  }
}
