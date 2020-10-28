function b = fstdeb_flipdims(a,n)
%FSTDEB_FLIPDIMS Debug function

[d1,d2]=size(a);
if mod(d1,n)~=0
  error('A wrong size');
end
d1=d1/n;
b=reshape(permute(reshape(a,n,d1,d2),[1 3 2]),n*d2,d1);
