function dg=fstdeb_diagmul(a,b,trs,varargin)
%FSTDEB_DIAGMUL Debug code

sa=fstdebug_getmat(a,1);
sb=fstdebug_getmat(b,1);
d=[];
if nargin>3
  d=varargin{1};
  if size(d,2)~=1
    d=d';
  end
end
if isempty(d)
  dg=sum(fstdebug_op(sa.mat,trs).*fstdebug_op(sb.mat,trs),1)';
else
  dg=sum(muldiag_old(d,fstdebug_op(sa.mat,trs)).*fstdebug_op(sb.mat, ...
						  trs),1)';
end
