function sa=fstdebug_getmat(a,varargin)
%FSTDEBUG_GETMAT Helper for debug code
%  SA = FSTDEBUG_GETMAT(A,{ignore=0})
%  A is matrix argument passed to FST_XXX function. Returns
%  structure SA, where SA.MAT is the matrix, SA.SPOS are the
%  starting positions [YS,XS], and SA.SCODE is the structure code
%  string. If A contains structure codes, the corr. pattern is
%  applied to SA.MAT. This is not done if IGNORE==1

sa=[];
ignore=0;
if nargin>1
  ignore=varargin{1};
end
if ~iscell(a)
  sa.mat=a;
  sa.spos=[1; 1];
  sa.scode='  ';
else
  if length(a)<2
    error('Error in A');
  end
  sa.spos=a{2};
  if length(sa.spos)~=4
    error('Error in A, SPOS entry');
  end
  ys=sa.spos(1); xs=sa.spos(2); m=sa.spos(3); n=sa.spos(4);
  sa.spos=[ys; xs];
  if ys<1 | xs<1 | m<0 | n<0 | ys+m-1>size(a{1},1) | xs+n-1> ...
	size(a{1},2)
    error('Error in A, SPOS entry');
  end
  sa.mat=a{1}(ys:(ys+m-1),xs:(xs+n-1));
  if length(a)>2
    sa.scode=a{3};
    if ~ischar(sa.scode) | length(sa.scode)~=2
      error('Error in A, SCODE entry');
    end
    suplo=sa.scode(1);
    if suplo~=' ' & suplo~='U' & suplo~='L'
      error('Error in A, SCODE(UPLO) entry');
    end
    sdiag=sa.scode(2);
    if sdiag~=' ' & sdiag~='U' & sdiag~='N'
      error('Error in A, SCODE(DIAG) entry');
    end
    if suplo==' ' & sdiag~=' '
      error('Error in A, SCODE entry');
    end
    if ~ignore
      if suplo~=' '
	if m~=n
	  error('Error in A, need square matrix');
	end
	if suplo=='L'
	  sa.mat=tril(sa.mat);
	else
	  sa.mat=triu(sa.mat);
	end
	if sdiag=='U'
	  sa.mat(1:(n+1):(n*n))=1;
	end
      end
    end
  else
    sa.scode='  ';
  end
end
