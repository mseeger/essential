function y = fstdeb_sumpos(x,ind,yn)
%FSTDEB_SUMPOS Debug function

y=zeros(yn,1);
n=length(ind); ind=reshape(ind,n,1);
pos=unique(ind)'; nu=length(pos);
mat=(ind(:,ones(nu,1))==pos(ones(n,1),:));
y(pos)=mat'*x;
