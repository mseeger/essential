function reta=fstdeb_muldiag(a,d,left)
%FSTDEB_MULDIAG Debug code

sa=fstdebug_getmat(a,1);
if left
  reta=fstdebug_writeback(a,muldiag_old(d,sa.mat),1);
else
  reta=fstdebug_writeback(a,muldiag_old(sa.mat,d),1);
end
