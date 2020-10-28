function mat=fstdebug_writeback(c,a,varargin)
%FSTDEBUG_WRITEBACK Helper for debug code
%  MAT = FSTDEBUG_WRITEBACK(C,A,{IGNORE=0})
%  A is a matrix. C is a matrix arg. which is a matrix or a cell
%  array. If C is a matrix, then MAT is set to A. If C is a cell
%  array, MAT is set equal to C{1}, then
%  the part corr. to C{2} is overwritten by A. If structure codes
%  are specified in C{3} and IGNORE==0, a structure pattern may be
%  applied for this copying.

ignore=0;
if nargin>2
  ignore=varargin{1};
end
if ~iscell(c)
  mat=a;
else
  mat=c{1};
  ys=c{2}(1); xs=c{2}(2); m=c{2}(3); n=c{2}(4);
  if length(c)<3 | c{3}(1)==' ' | ignore
    mat(ys:(ys+m-1),xs:(xs+n-1))=a;
  elseif c{3}(1)=='L'
    if c{3}(2)~='U'
      mat(ys:(ys+m-1),xs:(xs+n-1))=triu(mat(ys:(ys+m-1),xs:(xs+n- ...
						  1)),1)+tril(a);
    else
      mat(ys:(ys+m-1),xs:(xs+n-1))=triu(mat(ys:(ys+m-1),xs:(xs+n- ...
						  1)))+tril(a,-1);
    end
  else
    if c{3}(2)~='U'
      mat(ys:(ys+m-1),xs:(xs+n-1))=tril(mat(ys:(ys+m-1),xs:(xs+n- ...
						  1)),-1)+triu(a);
    else
      mat(ys:(ys+m-1),xs:(xs+n-1))=tril(mat(ys:(ys+m-1),xs:(xs+n- ...
						  1)))+triu(a,1);
    end
  end
end
