/* -------------------------------------------------------------------
 * FST_PERMUTE
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * Replace A by A(:,PERM)  (COL==1)
 * Replace A by A(PERM,:)  (COL==0)
 *
 * Structure codes ignored. Permutation is done in place, with
 * single col./row as temp. storage.
 *
 * Input:
 * - A:     Input/output matrix A
 * - PERM:  Permutation index
 * - COL:   Flag, s.a.
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

#define SETTARG(pos) do { \
  targ=((int) pvec[pos])-1; \
  if (targ<0 || targ>=n) \
    mexErrMsgTxt("PERM has invalid entries"); \
} while (0)

/* Main function FST_PERMUTE */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,m,n,lda,nxt,ione=1,start,numDone,targ;
  fst_matrix amat;
  const double* pvec;
  double* tempvec;
  bool* tickoff;
  bool col,emptyCycle;

  /* Read arguments */
  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  parseBLASMatrix(prhs[0],"A",&amat,-1,-1);
  col=(getScalInt(prhs[2],"COL")!=0);
  if (col) {
    m=amat.m; n=amat.n;
    lda=1; nxt=amat.stride;
  } else {
    m=amat.n; n=amat.m;
    lda=amat.stride; nxt=1;
  }
  if (getVecLen(prhs[1],"PERM")!=n)
    mexErrMsgTxt("PERM has wrong size");
  pvec=mxGetPr(prhs[1]);

  /* Temo. storage */
  tempvec=(double*) mxMalloc(m*sizeof(double));
  tickoff=(bool*) mxMalloc(n*sizeof(bool));
  for (i=0; i<n; i++) tickoff[i]=false;
  /* Main loop */
  start=0; numDone=0;
  for (;;) {
    /* Start new cycle with 'start' */
    BLASFUNC(dcopy) (&m,amat.buff+(start*nxt),&lda,tempvec,&ione);
    numDone++;
    SETTARG(start);
    emptyCycle=true;
    while (targ!=start) {
      BLASFUNC(dswap) (&m,amat.buff+(targ*nxt),&lda,tempvec,&ione);
      tickoff[targ]=true; numDone++; emptyCycle=false;
      SETTARG(targ);
      /* Sanity check (avoids infinite cycles) */
      if (numDone>n)
	mexErrMsgTxt("PERM is not a permutation");
    }
    if (!emptyCycle) {
      /* Close cycle */
      BLASFUNC(dcopy) (&m,tempvec,&ione,amat.buff+(start*nxt),&lda);
    }
    if (numDone==n) break; /* leave loop */
    /* Find next starting point */
    for (start++; start<n && tickoff[start]; start++);
    /* Sanity check */
    if (start==n)
      mexErrMsgTxt("PERM is not a permutation");
  }

  mxFree((void*) tempvec); mxFree((void*) tickoff);
}
