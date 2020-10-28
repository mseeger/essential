function y = gammacdf(x,alpha,omega)
%GAMMACDF Gamma cdf
%  GAMMACDF(X,ALPHA,OMEGA) computes the cdf of the Gamma distribution with
%  mean OMEGA and shape parameter ALPHA. The variance of the distribution is
%  2*OMEGA^2/ALPHA.
%  NOTE: If ALPHA is a natural number, it is the d.o.f. parameter, i.e. we
%  have a Chi-Square cdf.
y=gammainc((0.5*alpha/omega)*x,0.5*alpha);
