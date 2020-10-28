/* -------------------------------------------------------------------
 * FST_ACCUMULATE
 *
 * Specialized accumulation function. A is NN-by-NC matrix, IV an
 * NN index vector, B cell array of index vectors. For each i,
 * RET(i) is the sum of A(i,j), where j runs over B{IV(i)}.
 * If BTR is given, it is an additional index. In this case, RET(i)
 * is the sum of A(i,BTR(j)), where j runs over B{IV(i)}.
 *
 * Input:
 * - A:     Input matrix
 * - IV:    Index vector
 * - B:     Cell array of index vectors
 * - BTR:   Index vector, optional
 *
 * Return:
 * - RET:   Vector, s.a.
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include "mex.h"
#include "mex_helper.h"

char errMsg[200];

/* Main function FST_ */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,nn,nc,nb,ielem,sz,pos,s;
  double sum;
  fst_matrix amat;
  const double* iv,*bind,*aP,*btr=0;
  double* retP;
  const mxArray* bmat;
  int* ibtr=0;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=1)
    mexErrMsgTxt("Returns one argument");
  parseBLASMatrix(prhs[0],"A",&amat,-1,-1);
  nn=amat.m; nc=amat.n;
  if (getVecLen(prhs[1],"IV")!=nn)
    mexErrMsgTxt("IV has wrong size");
  iv=mxGetPr(prhs[1]);
  if (!mxIsCell(prhs[2]))
    mexErrMsgTxt("B must be cell array");
  if ((nb=mxGetM(prhs[2])*mxGetN(prhs[2]))==0)
    mexErrMsgTxt("B must not be empty");
  if (nrhs>3) {
    if (!mxIsDouble(prhs[3]))
      mexErrMsgTxt("Need vector in BTR");
    btr=mxGetPr(prhs[3]);
    sz=mxGetN(prhs[3])*mxGetM(prhs[3]);
    ibtr=(int*) mxMalloc(sz*sizeof(int));
    for (i=0; i<sz; i++) ibtr[i]=((int) btr[i])-1;
  }

  plhs[0]=mxCreateDoubleMatrix(nn,1,mxREAL);
  retP=mxGetPr(plhs[0]);
  aP=amat.buff;
  s=amat.stride;
  for (i=0; i<nn; i++,aP++) {
    ielem=((int) iv[i])-1;
    if (ielem<0 || ielem>=nb)
      mexErrMsgTxt("IV has wrong entries");
    bmat=mxGetCell(prhs[2],ielem);
    sz=mxGetN(bmat)*mxGetM(bmat);
    bind=mxGetPr(bmat);
    if (btr==0) {
      for (j=0,sum=0.0; j<sz; j++) {
	pos=((int) bind[j])-1;
	sum+=*(aP+(s*pos));
      }
    } else {
      for (j=0,sum=0.0; j<sz; j++) {
	pos=((int) bind[j])-1;
	sum+=*(aP+(s*ibtr[pos]));
      }
    }
    *(retP++)=sum;
  }
  if (ibtr!=0) mxFree((void*) ibtr);
}
