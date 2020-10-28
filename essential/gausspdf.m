function y=gausspdf(x)
%GAUSSPDF Density function of standard N(0,1) Gaussian
% Y=GAUSSPDF(X) returns standard normal p.d.f. at X

y=exp(-0.5*x.*x)/sqrt(2*pi);
