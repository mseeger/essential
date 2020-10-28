function retb=fstdeb_dtrmm(b,a,left,trs,varargin)
%FSTDEB_DTRMM Debug code

sb=fstdebug_getmat(b,1);
sa=fstdebug_getmat(a);
alpha=1;
if nargin>4
  alpha=varargin{1};
end
if sa.scode(1)==' '
  error('A must have UPLO');
end
if left
  retb=fstdebug_writeback(b,alpha*fstdebug_op(sa.mat,trs)*sb.mat, ...
			  1);
else
  retb=fstdebug_writeback(b,alpha*sb.mat*fstdebug_op(sa.mat,trs), ...
			  1);
end
