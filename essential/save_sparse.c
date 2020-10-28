/* -------------------------------------------------------------------
 * SAVE_SPARSE
 *
 * Stores sparse matrix A into text file with name FNAME. Here, 0
 * is stored as "0 ", requiring 2 bytes only. If VEC is given, it
 * is appended left to A, and its entries are stored as integers.
 *
 * Input:
 * - A:     Sparse input matrix
 * - FNAME: Name of file to write
 * - VEC:   S.a, optional
 * -------------------------------------------------------------------
 * Matlab MEX Function
 * Author: Matthias Seeger
 * ------------------------------------------------------------------- */

#define MATLAB_VER65

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "mex_helper.h"

char errMsg[200];

/* Main function SAVE_SPARSE */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,k,m,n;
  const double* amat,*vec=0;
  const int* icol,*irow;
  const char* fname;
  FILE* fid;
  int* pos,*rpos;

  /* Read arguments */
  if (nrhs<2)
    mexErrMsgTxt("Not enough input arguments");
  if (!mxIsSparse(prhs[0]))
    mexErrMsgTxt("A not sparse");
  m=mxGetM(prhs[0]); n=mxGetN(prhs[0]);
  amat=mxGetPr(prhs[0]); icol=mxGetJc(prhs[0]);
  irow=mxGetIr(prhs[0]);
  fname=getString(prhs[1],"FNAME");
  if (nrhs>2) {
    if (getVecLen(prhs[2],"VEC")!=m)
      mexErrMsgTxt("VEC has wrong size");
    vec=mxGetPr(prhs[2]);
  }

  if ((fid=fopen(fname,"w"))==0)
    mexErrMsgTxt("Cannot create file with name from FNAME");
  pos=(int*) mxMalloc(n*sizeof(int));
  rpos=(int*) mxMalloc(n*sizeof(int));
  for (j=0; j<n; j++) {
    pos[j]=k=icol[j];
    rpos[j]=(k<icol[j+1])?irow[k]:-1;
  }
  for (i=0; i<m; i++) {
    if (vec!=0)
      fprintf(fid,"%d ",(int) vec[i]);
    for (j=0; j<n; j++) {
      if (i==rpos[j]) {
	fprintf(fid,"%f ",amat[pos[j]]);
	k=++pos[j];
	rpos[j]=(k<icol[j+1])?irow[k]:-1;
      } else
	fprintf(fid,"0 ");
    }
    fprintf(fid,"\n");
  }
  fclose(fid);
  mxFree((void*) fname);
  mxFree((void*) pos); mxFree((void*) rpos);
}
