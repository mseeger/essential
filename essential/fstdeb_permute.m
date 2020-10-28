function reta=fstdeb_permute(a,perm,col)
%FSTDEB_PERMUTE Debug code

sa=fstdebug_getmat(a,1);
if col
  tempm(:,perm)=sa.mat;
  reta=fstdebug_writeback(a,tempm,1);
else
  tempm(perm,:)=sa.mat;
  reta=fstdebug_writeback(a,tempm,1);
end
