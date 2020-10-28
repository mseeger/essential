function mat = radialcf(s,t,c,w,varargin)
%RADIALCF RBF covariance function
%  RADIALCF(S,T,C,W,{VB}): S is N by D, T is M by D. C and W are positive
%  kernel parameters (below). Returns N by M matrix computed by applying
%  kernel to each pair of rows from S and T. The kernel is
%    K(x,y) = C*exp(-(0.5/D)*W*(x-y)'*(x-y)) + VB.
%  Here, VB==0 if not given, otherwise must be positive.

[n d]=size(s);
[m dummy]=size(t);
if dummy~=d
  error('Wrong dimension of inputs');
end
mat=exp(w*radialcf_precmat(s,t)+log(c));
if nargin>4
  vb=varargin{1};
  mat=mat+vb;
end
