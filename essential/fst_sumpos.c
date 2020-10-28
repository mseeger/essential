/* -------------------------------------------------------------------
 * FST_SUMPOS
 *
 * Y is obtained by starting from 0 of size YN, then adding X(j)
 * to Y(IND(j)) for all j.
 *
 * Input:
 * - X:   Input vector
 * - IND: Index
 * - YN:  S.a.
 * Return:
 * - Y:   Return vector
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include "mex.h"
#include "mex_helper.h"

char errMsg[200];

/* Main function FST_SUMPOS */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,yn,n,pos;
  const double* xvec,*ind;
  double* yvec;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=1)
    mexErrMsgTxt("Function returns 1 argument");
  if ((n=getVecLen(prhs[0],"X"))==0)
    mexErrMsgTxt("X is empty");
  xvec=mxGetPr(prhs[0]);
  if (getVecLen(prhs[1],"IND")!=n)
    mexErrMsgTxt("IND has wrong size");
  ind=mxGetPr(prhs[1]);
  yn=getScalInt(prhs[2],"YN");
  if (yn<1)
    mexErrMsgTxt("YN wrong");

  plhs[0]=mxCreateDoubleMatrix(yn,1,mxREAL);
  yvec=mxGetPr(plhs[0]);
  for (i=0; i<n; i++) {
    pos=((int) *(ind++))-1;
    if (pos<0 || pos>=yn)
      mexErrMsgTxt("IND or YN wrong");
    yvec[pos]+=(*(xvec++));
  }
}
