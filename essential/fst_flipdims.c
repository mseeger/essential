/* -------------------------------------------------------------------
 * FST_FLIPDIMS
 *
 * ATTENTION: We use the undocumented fact that the content of
 * matrices passed as arguments to a MEX function can be overwritten
 * like in a proper call-by-reference. This is not officially
 * supported and may not work in future Matlab versions!
 *
 * Operates on A, which is of size N*D1-by-D2. The flat buffer of
 * A is taken as N-by-D1-by-D2 tensor, and dimensions 2 and 3
 * are exchanged. The resulting matrix is reshaped to have D1
 * columns. In-place variant of
 *   RESHAPE(PERMUTE(RESHAPE(A,N,D1,D2),[1 3 2]),N*D2,D1),
 * but Matlab PERMUTE is not in-place.
 *
 * NOTE: This is not straightforward if D1~=D2. We use CACM Algorithm
 * 467 which does transposition in place:
 *   Brenner: Matrix Transposition in Place, CACM 16(11), 1973
 *   Available: portal.acm.org
 * See also matrix/TransInPlace.h of STATSIM.
 *
 * Input:
 * - A:       Matrix, overwritten by result
 * - N:       Size of first dimension
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

/* Code ripped off from matrix/TransInPlace.h of STATSIM. See comments
   there. */

int factor(int n,int* ifact,int* ipower,int* nexp)
{
  int ip=-1,ifcur=0,npart=n,idiv=2;
  int iquot;

  for (;;) {
    iquot=npart/idiv;
    if ((npart-iquot*idiv)==0) {
      if (idiv>ifcur) {
	ip++;
	ifact[ip]=ipower[ip]=idiv; nexp[ip]=1;
	ifcur=idiv;
      } else {
	ipower[ip]*=idiv; nexp[ip]++;
      }
      npart=iquot;
    } else {
      if (iquot<=idiv) break;
      if (idiv==2)
	idiv=3;
      else
	idiv+=2;
    }
  }
  if (npart>1) {
    if (npart>ifcur) {
      ip++;
      ifact[ip]=ipower[ip]=npart; nexp[ip]=1;
    } else {
      ipower[ip]*=npart; nexp[ip]++;
    }
  }
  return (ip+1);
}

/* Main function FST_FLIPDIMS */

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j,ip,nn,n1,n2,ione=1;
  int n,n12,m;
  double temp;
  double* amat;
  int num,off;
  double* a1P,*a2P;
  int nwork,npower;
  bool* moved;
  int ifact[8],ipower[8],nexp[8],iexp[8];
  int idiv,ncount,istart,mmist,isoid,itest,itemp;
  int ia1,ia2,mmia1,mmia2;
  bool cont;
  double* atemp,*btemp;

  /* Read arguments */
  if (nrhs<2)
    mexErrMsgTxt("Not enough input arguments");
  if ((nn=getScalInt(prhs[1],"N"))<=0)
    mexErrMsgTxt("Wrong argument N");
  checkMatrix(prhs[0],"A",-1,-1);
  if ((n1=mxGetM(prhs[0]))%nn!=0)
    mexErrMsgTxt("A has wrong size");
  amat=mxGetPr(prhs[0]);
  n1/=nn; n2=mxGetN(prhs[0]);
  n12=n1*n2; m=n12-1; n=n1;

  if (n1>1 && n2>1) {
    /* Transpose in place */
    if (n1==n2) {
      /* Square matrix: simple special case */
      a1P=amat;
      for (i=1,off=n; i<=n; i++,off+=(n+1)) {
	a1P+=nn*i; a2P=amat+nn*off;
	num=n-i;
	for (j=0; j<num; j++) {
	  if (nn>1)
	    BLASFUNC(dswap) (&nn,a1P,&ione,a2P,&ione);
	  else {
	    temp=*a1P; *a1P=*a1P; *a2P=temp;
	  }
	  a1P+=nn; a2P+=nn*n;
	}
      }
    } else {
      /* General rectangular matrix */
      nwork=(n1+n2)/2;
      moved=(bool*) mxMalloc(nwork*sizeof(bool));
      npower=factor(m,ifact,ipower,nexp);
      for (ip=0; ip<npower; ip++) iexp[ip]=0;
      /* Generate every divisor of m less than m/2 */
      idiv=1;
      atemp=(double*) mxMalloc(nn*sizeof(double));
      btemp=(double*) mxMalloc(nn*sizeof(double));
      while (idiv<m/2) {
	/* The number of elements whose index is div. by idiv and by no other
	   divisor of m is the Euler totient function, phi(m/idiv) */
	ncount=m/idiv;
	for (ip=0; ip<npower; ip++)
	  if (iexp[ip]!=nexp[ip])
	    ncount=(ncount/ifact[ip])*(ifact[ip]-1);
	for (i=0; i<nwork; i++) moved[i]=false;
	/* The starting point of a subcycle is div. only by idiv and must
	   not appear in any other subcycle */
	for (istart=idiv; ncount>0; istart+=idiv) {
	  mmist=m-istart;
	  if (istart!=idiv) {
	    if (istart<=nwork && moved[istart-1]) continue;
	    isoid=istart/idiv;
	    for (ip=0,cont=false; ip<npower; ip++)
	      if (iexp[ip]!=nexp[ip] && (isoid%ifact[ip])==0) {
		cont=true; break;
	      }
	    if (cont) continue;
	    if (istart>nwork) {
	      itest=istart; cont=false;
	      do {
		itemp=itest/n2;
		itest=itemp+(itest-itemp*n2)*n; /* succ. in permut. */
		if (itest<istart || itest>mmist) {
		  cont=true; break;
		}
	      } while (itest>istart && itest<mmist);
	      if (cont) continue;
	    }
	  }
	  /* atemp=a[istart]; btemp=a[mmist]; */
	  if (nn>1) {
	    BLASFUNC(dcopy) (&nn,amat+nn*istart,&ione,atemp,&ione);
	    BLASFUNC(dcopy) (&nn,amat+nn*mmist,&ione,btemp,&ione);
	  } else {
	    *atemp=amat[istart]; *btemp=amat[mmist];
	  }
	  for (ia1=istart; ; ) {
	    itemp=ia1/n2;
	    ia2=itemp+(ia1-itemp*n2)*n; /* succ. in permut. */
	    mmia1=m-ia1; mmia2=m-ia2;
	    if (ia1<=nwork) moved[ia1-1]=true;
	    if (mmia1<=nwork) moved[mmia1-1]=true;
	    ncount-=2;
	    /* Move two elements, the second from the negative subcycle.
	       Check first for subcycle closure */
	    if (ia2==istart) {
	      /* a[ia1]=atemp; a[mmia1]=btemp; */
	      if (nn>1) {
		BLASFUNC(dcopy) (&nn,atemp,&ione,amat+nn*ia1,&ione);
		BLASFUNC(dcopy) (&nn,btemp,&ione,amat+nn*mmia1,&ione);
	      } else {
		amat[ia1]=*atemp; amat[mmia1]=*btemp;
	      }
	      break;
	    } else if (mmia2==istart) {
	      /* a[ia1]=btemp; a[mmia1]=atemp; */
	      if (nn>1) {
		BLASFUNC(dcopy) (&nn,btemp,&ione,amat+nn*ia1,&ione);
		BLASFUNC(dcopy) (&nn,atemp,&ione,amat+nn*mmia1,&ione);
	      } else {
		amat[ia1]=*btemp; amat[mmia1]=*atemp;
	      }
	      break;
	    } else {
	      /* a[ia1]=a[ia2]; a[mmia1]=a[mmia2]; */
	      if (nn>1) {
		BLASFUNC(dcopy) (&nn,amat+nn*ia2,&ione,amat+nn*ia1,&ione);
		BLASFUNC(dcopy) (&nn,amat+nn*mmia2,&ione,amat+nn*mmia1,&ione);
	      } else {
		amat[ia1]=amat[ia2]; amat[mmia1]=amat[mmia2];
	      }
	      ia1=ia2;
	    }
	  }
	}
	for (ip=0,cont=false; ip<npower; ip++) {
	  if (iexp[ip]!=nexp[ip]) {
	    iexp[ip]++;
	    idiv*=ifact[ip];
	    cont=true; break;
	  } else {
	    iexp[ip]=0;
	    idiv/=ipower[ip];
	  }
	}
	if (!cont) break; /* leave loop */
      }
      mxFree((void*) moved);
      mxFree((void*) atemp); mxFree((void*) btemp);
    }
  }
  /* Reshape */
  if (n1<=n2) {
    mxSetN((mxArray*) prhs[0],n1);
    mxSetM((mxArray*) prhs[0],nn*n2);
  } else {
    mxSetM((mxArray*) prhs[0],nn*n2);
    mxSetN((mxArray*) prhs[0],n1);
  }
}
