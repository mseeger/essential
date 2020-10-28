function y = gammadens(x,shape,iscale)
%GAMMADENS Gamma density
% Y=GAMMADENS(X,SHAPE,ISCALE) computes density of the Gamma distribution with
%  shape SHAPE and inverse scale ISCALE. The density is
%    f(x) = [(b/2)^{a/2}/Gamma(a/2)]*x^{a/2-1}*exp(-(b/2)*x),
%  where a is the shape, b the inverse scale parameter. The mean is a/b, the
%  variance is 2*a/b^2.
%  Note: For a=nu,b=1, we have the Chi-Square density with nu d.o.f.

p1=shape/2;
p2=iscale/2;
offset=-gammaln(p1)+p1*log(p2);
p1=p1-1;
y=exp(offset+p1*log(x)-p2*x).*(x>0);
