/* -------------------------------------------------------------------
 * FST_GETSPARSE
 *
 * Returns the Pr, Ir, and Jc arrays for a sparse matrix A.
 *
 * Input:
 * - A:  Sparse input matrix A
 * Return:
 * - PR: Pr array
 * - IR: Ir array (int32)
 * - JC: Jc array (int32)
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

/* Main function FST_GETSPARSE */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int sz,n;
  const int* jc;
  int dims[2];

  /* Sanity checks */
  if (sizeof(double)!=8 || sizeof(int)!=4)
    mexErrMsgTxt("Need byte size 4 for int, 8 for double");
  /* Read arguments */
  if (nrhs<1)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=3)
    mexErrMsgTxt("Returns 3 arguments");
  if (!mxIsSparse(prhs[0]))
    mexErrMsgTxt("A not sparse matrix");
  n=mxGetN(prhs[0]);
  jc=mxGetJc(prhs[0]);
  sz=jc[n];
  plhs[0]=mxCreateDoubleMatrix(sz,1,mxREAL);
  memmove((char*) mxGetPr(plhs[0]),(char*) mxGetPr(prhs[0]),8*sz);
  dims[0]=sz; dims[1]=1;
  plhs[1]=mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
  memmove((char*) mxGetData(plhs[1]),(char*) mxGetIr(prhs[0]),4*sz);
  dims[0]=n+1;
  plhs[2]=mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
  memmove((char*) mxGetData(plhs[2]),(char*) jc,4*(n+1));
}
