function sa=fstdebug_sampmat(m,n,varargin)
%FSTDEBUG_SAMPMAT Samples random matrix
%  SA = FSTDEBUG_SAMPMAT(M,N,{PARTONLY=0})
%
%  With prob. 2/3, the matrix is part of a larger buffer matrix,
%  whose number of rows/cols can be up to twice the matrix size.
%  If PARTONLY==1, the matrix is always a part.

partonly=0;
if nargin>2
  partonly=varargin{1};
end
if ~partonly & rand<(1/3)
  sa=2*rand(m,n)-1;
else
  bm=m+1+floor(rand*m);
  bn=n+1+floor(rand*n);
  ys=1+floor(rand*(bm-m)); xs=1+floor(rand*(bn-n));
  sa{1}=2*rand(bm,bn)-1;
  sa{2}=[ys; xs; m; n];
  sa{3}='  ';
end
