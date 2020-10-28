/* -------------------------------------------------------------------
 * FST_RESHAPE
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * In-place version of RESHAPE. A is a matrix whose buffer has ELA
 * elements. It is reshaped to size M-by-N, provided that
 * ELA >= M*N. If ELA > M*N, the elements beyond M*N are lost.
 *
 * Input:
 * - A:   Input/output matrix
 * - M:   New number rows
 * - N:   New number cols
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include "mex.h"
#include "mex_helper.h"

char errMsg[200];

/* Main function FST_RESHAPE */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int m,n;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  checkMatrix(prhs[0],"A",-1,-1);
  m=getScalInt(prhs[1],"M");
  n=getScalInt(prhs[2],"N");
  if (m<1 || n<1)
    mexErrMsgTxt("Wrong M, N");
  if (mxGetNumberOfElements(prhs[0])<m*n)
    mexErrMsgTxt("A is too small");
  if (m<=mxGetM(prhs[0])) {
    mxSetM((mxArray*) prhs[0],m);
    mxSetN((mxArray*) prhs[0],n);
  } else {
    mxSetN((mxArray*) prhs[0],n);
    mxSetM((mxArray*) prhs[0],m);
  }
}
