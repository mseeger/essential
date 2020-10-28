function retc=fstdeb_dspamm(c,a,b,tra,varargin)
%FSTDEB_DSPAMM Debug function

sc=fstdebug_getmat(c,1);
sb=fstdebug_getmat(b,1);
if ~issparse(a)
  error('A must be sparse');
end
alpha=1; beta=0;
if nargin>4
  alpha=varargin{1};
  if nargin>5
    beta=varargin{2};
  end
end
if tra
  retc=fstdebug_writeback(c,beta*sc.mat+full((alpha*a')*sb.mat),1);
else
  retc=fstdebug_writeback(c,beta*sc.mat+full((alpha*a)*sb.mat),1);
end
