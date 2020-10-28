function y=logsumexp(x)
%LOGSUMEXP Sums exponentiated comps. and takes log, in a stable way
%  Y=LOGSUMEXP(X). Exponentiates each column of X, sums the components and
%  computes the log. This is done in a stable way.
%  NOTE: The implementation temp. creates a matrix of size of X.

mx=max(x);
y=mx+log(sum(exp(x-mx(ones(size(x,1),1),:)),1));
