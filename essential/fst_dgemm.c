/* -------------------------------------------------------------------
 * FST_DGEMM
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*op(A)*op(B) + BETA*C, op(X) == X or X'
 *
 * Input:
 * - C:     Input/output matrix. Must have correct size, even if
 *          BETA==0!
 * - A:     Input matrix A
 * - TRA:   Use A' for A?
 * - B:     Input matrix B
 * - TRB:   Use B' for B?
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
  fst_matrix amat,bmat,cmat;
  char tra[2],trb[2];

  /* Read arguments */
  if (nrhs<5)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  tra[1]=0;
  if (getScalInt(prhs[2],"TRA")!=0) {
    tra[0]='T';
    parseBLASMatrix(prhs[1],"A",&amat,-1,cmat.m);
    k=amat.m;
  } else {
    tra[0]='N';
    parseBLASMatrix(prhs[1],"A",&amat,cmat.m,-1);
    k=amat.n;
  }
  trb[1]=0;
  if (getScalInt(prhs[4],"TRB")!=0) {
    trb[0]='T';
    parseBLASMatrix(prhs[3],"B",&bmat,cmat.n,k);
  } else {
    trb[0]='N';
    parseBLASMatrix(prhs[3],"B",&bmat,k,cmat.n);
  }
  alpha=1.0; beta=0.0;
  if (nrhs>5) {
    alpha=getScalar(prhs[5],"ALPHA");
    if (nrhs>6)
      beta=getScalar(prhs[6],"BETA");
  }

  BLASFUNC(dgemm) (tra,trb,&cmat.m,&cmat.n,&k,&alpha,amat.buff,&amat.stride,
		   bmat.buff,&bmat.stride,&beta,cmat.buff,&cmat.stride);
}
