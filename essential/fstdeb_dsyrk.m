function retc=fstdeb_dsyrk(c,a,trs,varargin)
%FSTDEB_DSYRK Debug for FST_DSYRK

sc=fstdebug_getmat(c);
sa=fstdebug_getmat(a,1);
alpha=1; beta=0;
if nargin>3
  alpha=varargin{1};
  if nargin>4
    beta=varargin{2};
  end
end
if sc.scode(1)==' '
  error('C must have UPLO');
end
retc=fstdebug_writeback(c,alpha*fstdebug_op(sa.mat,trs)* ...
			fstdebug_op(sa.mat',trs)+beta*sc.mat);
