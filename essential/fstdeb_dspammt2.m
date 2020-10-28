function retd=fstdeb_dspammt2(d,a,b,c,varargin)
%FSTDEB_DSPAMMT2 Debug function

sd=fstdebug_getmat(d,1);
sc=fstdebug_getmat(c,1);
if ~issparse(a)
  error('A must be sparse');
end
if ~issparse(b)
  error('B must be sparse');
end
alpha=1; beta=0;
if nargin>4
  alpha=varargin{1};
  if nargin>5
    beta=varargin{2};
  end
end
retd=fstdebug_writeback(d,beta*sd.mat+full((alpha*a)*full(b'*sc.mat)),1);
