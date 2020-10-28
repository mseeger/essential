%FST_DTRSM Operator backsubstitution with triangular matrix
%  FST_DTRSM(B,A,LEFT,TRS,{ALPHA=1})
%
%  B = ALPHA*op(A^-1)*B (LEFT==1), op(X) == X or X', dep. on TRS
%  B = ALPHA*B*op(A^-1) (LEFT==0)
% 
%  Here, A must be triangular (UPLO str. code required, DIAG code
%  optional) and nonsingular. If B has str. codes, they are
%  ignored.
