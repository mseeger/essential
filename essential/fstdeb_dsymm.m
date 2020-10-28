function retc=fstdeb_dsymm(c,a,b,varargin)
%FSTDEB_DSYMM Debug for FST_DSYMM

sc=fstdebug_getmat(c,1);
sa=fstdebug_getmat(a);
sb=fstdebug_getmat(b);
alpha=1; beta=0;
if nargin>3
  alpha=varargin{1};
  if nargin>4
    beta=varargin{2};
  end
end
if sa.scode(1)=='L'
  sa.mat=makesymm(sa.mat,1);
elseif sa.scode(1)=='U'
  sa.mat=makesymm(sa.mat,0);
elseif sb.scode(1)=='L'
  sb.mat=makesymm(sb.mat,1);
elseif sb.scode(1)=='U'
  sb.mat=makesymm(sb.mat,0);
else
  error('One of A, B must have UPLO');
end
retc=fstdebug_writeback(c,alpha*sa.mat*sb.mat+beta*sc.mat,1);
