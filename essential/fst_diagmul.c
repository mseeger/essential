/* -------------------------------------------------------------------
 * FST_DIAGMUL
 *
 * DG = DIAG(op(A)'*DIAG(D)*op(B)), op(X) == X or X', dep. on TRS
 *
 * A, B must have the same size. D is all ones by def. Structure
 * codes ignored.
 *
 * Input:
 * - A:     Input matrix A
 * - B:     Input matrix B
 * - TRS:   Flag, s.a.
 * - D:     S.a. Def.: all ones
 *
 * Return:
 * - DG:    Return vector
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

/* Predeclarations */

/* Main function FST_DIAGMUL */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,ione=1;
  double sum;
  fst_matrix amat,bmat;
  double* dgvec;
  const double* dvec=0,*aP,*bP,*dP;
  bool trs;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  if (nlhs!=1)
    mexErrMsgTxt("Has one return argument");
  parseBLASMatrix(prhs[0],"A",&amat,-1,-1);
  parseBLASMatrix(prhs[1],"B",&bmat,amat.m,amat.n);
  trs=(getScalInt(prhs[2],"TRS")!=0);
  if (nrhs>3) {
    i=getVecLen(prhs[3],"D");
    if ((trs && i!=amat.n) || (!trs && i!=amat.m))
      mexErrMsgTxt("D has wrong size");
    dvec=mxGetPr(prhs[3]);
  }

  plhs[0]=mxCreateDoubleMatrix(trs?amat.m:amat.n,1,mxREAL);
  dgvec=mxGetPr(plhs[0]);
  if (dvec==0) {
    if (!trs) {
      for (i=0; i<amat.n; i++)
	dgvec[i]=BLASFUNC(ddot) (&amat.m,amat.buff+(i*amat.stride),&ione,
				 bmat.buff+(i*bmat.stride),&ione);
    } else {
      for (i=0; i<amat.m; i++)
	dgvec[i]=BLASFUNC(ddot) (&amat.n,amat.buff+i,&amat.stride,
				 bmat.buff+i,&bmat.stride);
    }
  } else {
    if (!trs) {
      for (i=0; i<amat.n; i++) {
	aP=amat.buff+(i*amat.stride); bP=bmat.buff+(i*bmat.stride);
	dP=dvec;
	for (j=0,sum=0.0; j<amat.m; j++)
	  sum+=(*(dP++))*(*(aP++))*(*(bP++));
	dgvec[i]=sum;
      }
    } else {
      for (i=0; i<amat.m; i++) {
	aP=amat.buff+i; bP=bmat.buff+i; dP=dvec;
	for (j=0,sum=0.0; j<amat.n; j++,aP+=amat.stride,bP+=bmat.stride)
	  sum+=(*(dP++))*(*aP)*(*bP);
	dgvec[i]=sum;
      }
    }
  }
}
