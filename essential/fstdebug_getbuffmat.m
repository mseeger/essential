function a=fstdebug_getbuffmat(sa)

if ~iscell(sa)
  a=sa;
else
  a=sa{1};
end
