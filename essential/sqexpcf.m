function mat = sqexpcf(s,t,c,w,varargin)
%SQEXPCF Squared-exponential covariance function.
% SQEXPCF(S,T,C,W,{VB}): S is N by D, T is M by D. C is pos. scalar, W
%   pos. D-vector. If given, VB is pos. scalar. Returns N by M matrix
%   computed by applying kernel to each pair of rows from S and T. The
%   kernel is
%     K(x,y) = C*exp(-(0.5/D)*sum_i(w_i*(s_i-t_i)^2)) + VB,
%   where w_i are the components of W. The def. for VB is 0.

if nargin>4,
  vb=varargin{1};
else
  vb=0;
end
[n d]=size(s);
[m dummy]=size(t);
if dummy~=d
  error('Wrong dimension of inputs');
end
if length(w)~=d
  error('Wrong dimension of W');
end
w=reshape(w,d,1);

if issparse(s) | issparse(t)
  wt=muldiag(t,w);
  fact=-0.5/d;
  tvec=fact*sum(t.*wt,2)';
  svec=fact*sum(muldiag(s,w).*s,2);
  mat=full(((1/d)*s)*wt')+svec(:,ones(m,1))+tvec(ones(n,1),:);
else
  wt=t(:,:);
  fst_muldiag(wt,w,0);
  fact=-0.5/d;
  tvec=fact*fst_diagmul(t,wt,1);
  svec=fact*fst_diagmul(s,s,1,w);
  mat=zeros(n,m);
  fst_dgemm(mat,s,0,wt,1,1/d);
  fst_addvec(mat,svec,1);
  fst_addvec(mat,tvec,0);
end
mat=exp(mat+log(c));
if vb~=0
  mat=mat+vb;
end
