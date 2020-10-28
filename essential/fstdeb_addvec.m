function retc=fstdeb_addvec(c,x,col,varargin)
%FSTDEB_ADDVEC Debug function

sc=fstdebug_getmat(c,1);
sz=length(x); x=reshape(x,sz,1);
if sz~=size(sc.mat,1+(~col))
  error('X has wrong size');
end
alpha=1;
if nargin>3
  alpha=varargin{1};
end
tempx=alpha*x;
if col
  retc=fstdebug_writeback(c,sc.mat+tempx(:,ones(size(sc.mat,2),1)), ...
			  1);
else
  tempx=tempx';
  retc=fstdebug_writeback(c,sc.mat+tempx(ones(size(sc.mat,1),1),:), ...
			  1);
end
