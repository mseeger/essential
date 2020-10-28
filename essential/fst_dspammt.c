/* -------------------------------------------------------------------
 * FST_DSPAMMT
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*A*A'*B + BETA*C
 *
 * Here, A is a sparse matrix, while C, B are dense matrices.
 * Structure codes are ignored. Intermed. storage of size COLS(B)
 * is used, but A'*B does not have to be stored.
 *
 * Input:
 * - C:     Dense output matrix. Must have correct size!
 * - A:     Sparse input matrix A
 * - B:     Dense input matrix B
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

/* Main function FST_DSPAMM */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,ione=1,m,n,p,ck;
  double alpha,beta,temp;
  fst_matrix bmat,cmat;
  const double* amat;
  double* tvec;
  const int* icol,*irow;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  p=cmat.n; m=cmat.m;
  parseBLASMatrix(prhs[2],"B",&bmat,m,p);
  if (!mxIsSparse(prhs[1]) || mxGetM(prhs[1])!=m)
    mexErrMsgTxt("A not sparse or wrong size");
  n=mxGetN(prhs[1]);
  amat=mxGetPr(prhs[1]); icol=mxGetJc(prhs[1]);
  irow=mxGetIr(prhs[1]);
  alpha=1.0; beta=0.0;
  if (nrhs>3) {
    alpha=getScalar(prhs[3],"ALPHA");
    if (nrhs>4)
      beta=getScalar(prhs[4],"BETA");
  }

  tvec=(double*) mxMalloc(p*sizeof(double));
  for (j=0; j<p; j++)
    BLASFUNC(dscal) (&m,&beta,cmat.buff+(j*cmat.stride),&ione);
  for (j=0; j<n; j++) {
    temp=0.0;
    BLASFUNC(dscal) (&p,&temp,tvec,&ione);
    ck=icol[j+1];
    for (i=icol[j]; i<ck; i++)
      BLASFUNC(daxpy) (&p,amat+i,bmat.buff+irow[i],&bmat.stride,tvec,&ione);
    for (i=icol[j]; i<ck; i++) {
      temp=alpha*amat[i];
      BLASFUNC(daxpy) (&p,&temp,tvec,&ione,cmat.buff+irow[i],&cmat.stride);
    }
  }
  mxFree((void*) tvec);
}
