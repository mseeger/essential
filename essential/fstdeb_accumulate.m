function ret=fstdeb_accumulate(a,iv,b)
%FSTDEB_ACCUMULATE Debug function

sa=fstdebug_getmat(a,1);
nn=size(sa.mat,1);
ret=zeros(nn,1);
for i=1:nn
  ret(i)=sum(sa.mat(i,b{iv(i)}));
end
