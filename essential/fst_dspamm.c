/* -------------------------------------------------------------------
 * FST_DSPAMM
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * C = ALPHA*op(A)*B + BETA*C, op(A)==A / A', dep. on TRA
 *
 * Here, A is a sparse matrix, while C, B are dense matrices.
 * Structure codes are ignored.
 *
 * Input:
 * - C:     Dense output matrix. Must have correct size!
 * - A:     Sparse input matrix A
 * - B:     Dense input matrix B
 * - TRA:   Flag, s.a.
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
  const double* amat,*bP;
  double* cP;
  const int* icol,*irow;
  bool tra;

  /* Read arguments */
  if (nrhs<4)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"C",&cmat,-1,-1);
  p=cmat.n;
  parseBLASMatrix(prhs[2],"B",&bmat,-1,p);
  tra=(getScalInt(prhs[3],"TRA")!=0);
  if (tra) {
    m=bmat.m; n=cmat.m;
  } else {
    m=cmat.m; n=bmat.m;
  }
  if (!mxIsSparse(prhs[1]) || mxGetM(prhs[1])!=m || mxGetN(prhs[1])!=n)
    mexErrMsgTxt("A not sparse or wrong size");
  amat=mxGetPr(prhs[1]); icol=mxGetJc(prhs[1]);
  irow=mxGetIr(prhs[1]);
  alpha=1.0; beta=0.0;
  if (nrhs>4) {
    alpha=getScalar(prhs[4],"ALPHA");
    if (nrhs>5)
      beta=getScalar(prhs[5],"BETA");
  }

  for (j=0; j<p; j++)
    BLASFUNC(dscal) (&cmat.m,&beta,cmat.buff+(j*cmat.stride),&ione);
  if (!tra) {
    bP=bmat.buff;
    for (j=0; j<n; j++,bP++) {
      for (i=icol[j],ck=icol[j+1]; i<ck; i++) {
	temp=alpha*amat[i];
	BLASFUNC(daxpy) (&p,&temp,bP,&bmat.stride,cmat.buff+irow[i],
			 &cmat.stride);
      }
    }
  } else {
    cP=cmat.buff;
    for (j=0; j<n; j++,cP++) {
      for (i=icol[j],ck=icol[j+1]; i<ck; i++) {
	temp=alpha*amat[i];
	BLASFUNC(daxpy) (&p,&temp,bmat.buff+irow[i],&bmat.stride,cP,
			 &cmat.stride);
      }
    }
  }
}
