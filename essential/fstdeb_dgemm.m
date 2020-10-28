function retc=fstdeb_dgemm(c,a,tra,b,trb,varargin)
%FSTDEB_DGEMM Debug for FST_DGEMM

sc=fstdebug_getmat(c,1);
sa=fstdebug_getmat(a,1);
sb=fstdebug_getmat(b,1);
alpha=1; beta=0;
if nargin>5
  alpha=varargin{1};
  if nargin>6
    beta=varargin{2};
  end
end
retc=fstdebug_writeback(c,alpha*fstdebug_op(sa.mat,tra)* ...
			fstdebug_op(sb.mat,trb)+beta*sc.mat,1);
