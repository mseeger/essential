function c=muldiag(a,b)
%MULDIAG Multiplies matrix with diagonal matrix
%  C=MULDIAG(A,B) computes A*B where at least one of A,B must be a
%  diagonal matrix, repr. as a column vector.
%  NOTE: If one of A, B is a scalar, it has priority as diagonal
%  matrix: if A col. vec., B scalar, B is treated as diag. matrix.

[m,n]=size(a);
[p,q]=size(b);
if n==1
  if p==1 & q==1
    % Special case: B scalar
    c=b*a;
    return
  elseif m~=p
    error('Wrong dimensions!');
  end
  if ~issparse(b)
    c=b(:,:);
    fst_muldiag(c,full(a),1);
  else
    c=spdiags(full(a),0,m,m)*b;
  end
elseif q==1
  if n~=p
    error('Wrong dimensions!');
  end
  if ~issparse(a)
    c=a(:,:);
    fst_muldiag(c,full(b),0);
  else
    c=a*spdiags(full(b),0,n,n);
  end
else
  error('One of A,B must be a diagonal matrix (column vector)!');
end
