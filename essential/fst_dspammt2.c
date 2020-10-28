/* -------------------------------------------------------------------
 * FST_DSPAMMT2
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * D = ALPHA*A*B'*C + BETA*D
 *
 * Here, A, B are sparse matrices with the same number of cols, while
 * D, C are dense matrices. Structure codes are ignored. Intermed.
 * storage of size COLS(C) is used, but B'*C does not have to be
 * stored.
 *
 * Input:
 * - D:     Dense output matrix. Must have correct size!
 * - A:     Sparse input matrix A
 * - B:     Sparse input matrix B
 * - C:     Dense input matrix C
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

/* Main function FST_DSPAMMT2 */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,ione=1,n,p,ck;
  double alpha,beta,temp;
  fst_matrix cmat,dmat;
  const double* amat,*bmat;
  double* tvec;
  const int* iacol,*iarow,*ibcol,*ibrow;

  /* Read arguments */
  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"D",&dmat,-1,-1);
  p=dmat.n;
  parseBLASMatrix(prhs[3],"C",&cmat,-1,p);
  if (!mxIsSparse(prhs[1]) || mxGetM(prhs[1])!=dmat.m)
    mexErrMsgTxt("A not sparse or wrong size");
  n=mxGetN(prhs[1]);
  amat=mxGetPr(prhs[1]); iacol=mxGetJc(prhs[1]);
  iarow=mxGetIr(prhs[1]);
  if (!mxIsSparse(prhs[2]) || mxGetM(prhs[2])!=cmat.m ||
      mxGetN(prhs[2])!=n)
    mexErrMsgTxt("B not sparse or wrong size");
  bmat=mxGetPr(prhs[2]); ibcol=mxGetJc(prhs[2]);
  ibrow=mxGetIr(prhs[2]);
  alpha=1.0; beta=0.0;
  if (nrhs>4) {
    alpha=getScalar(prhs[4],"ALPHA");
    if (nrhs>5)
      beta=getScalar(prhs[5],"BETA");
  }

  tvec=(double*) mxMalloc(p*sizeof(double));
  for (j=0; j<p; j++)
    BLASFUNC(dscal) (&dmat.m,&beta,dmat.buff+(j*dmat.stride),&ione);
  for (j=0; j<n; j++) {
    temp=0.0;
    BLASFUNC(dscal) (&p,&temp,tvec,&ione);
    for (i=ibcol[j],ck=ibcol[j+1]; i<ck; i++)
      BLASFUNC(daxpy) (&p,bmat+i,cmat.buff+ibrow[i],&cmat.stride,tvec,&ione);
    for (i=iacol[j],ck=iacol[j+1]; i<ck; i++) {
      temp=alpha*amat[i];
      BLASFUNC(daxpy) (&p,&temp,tvec,&ione,dmat.buff+iarow[i],&dmat.stride);
    }
  }
  mxFree((void*) tvec);
}
