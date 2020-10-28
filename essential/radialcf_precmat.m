function mat = radialcf_precmat(s,t)
%RADIALCF_PRECMAT Precomputation matrix for RBF covariance function
%  RADIALCF_PRECMAT(S,T): S is N by D, T is M by D. See RADIALCF
%  for RBF kernel. Here, compute MAT s.t. C*exp(W*MAT) + VB would
%  be the kernel matrix.

[n d]=size(s);
[m dummy]=size(t);
if dummy~=d
  error('Wrong dimension of inputs');
end
fact=-1/(2*d);
if issparse(s)
  if ~issparse(t)
    error('S, T must either be both sparse or none');
  end
  svec=fact*full(sum(s.*s,2)); tvec=fact*full(sum(t.*t,2))';
  mat=full(((1/d)*s)*t')+svec(:,ones(m,1))+tvec(ones(n,1),:);
else
  svec=fact*fst_diagmul(s,s,1);
  tvec=fact*fst_diagmul(t,t,1);
  mat=zeros(n,m);
  fst_dgemm(mat,s,0,t,1,1/d);
  fst_addvec(mat,svec,1);
  fst_addvec(mat,tvec,0);
end
