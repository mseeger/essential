%FST_DIAGMUL Diagonal of matrix-matrix multiplication
%  DG = FST_DIAGMUL(A,B,TRS,{D})
%
%  DG = DIAG(op(A)'*DIAG(D)*op(B)), op(X) == X or X', dep. on TRS
%
%  A, B must have the same size. If D is not given, it is all ones.
%  Structure codes are ignored.
