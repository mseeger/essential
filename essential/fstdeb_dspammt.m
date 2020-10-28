function retc=fstdeb_dspammt(c,a,b,varargin)
%FSTDEB_DSPAMMT Debug function

sc=fstdebug_getmat(c,1);
sb=fstdebug_getmat(b,1);
if ~issparse(a)
  error('A must be sparse');
end
alpha=1; beta=0;
if nargin>3
  alpha=varargin{1};
  if nargin>4
    beta=varargin{2};
  end
end
retc=fstdebug_writeback(c,beta*sc.mat+full((alpha*a)*full(a'*sb.mat)),1);
