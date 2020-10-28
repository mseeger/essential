function a = makesymm(b,varargin)
%MAKESYMM Symmetrizes matrix from lower/upper triangle
%  A = MAKESYMM(B,{LOW=true}) returns in A the symmetric matrix which has
%  the same lower triangle and diagonal as B. If LOW==false, the upper
%  triangle is used instead.

if nargin<2 | varargin{1}
  a=tril(b)+tril(b,-1)';
else
  a=triu(b)+triu(b,1)';
end
