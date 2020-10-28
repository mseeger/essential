%FST_DTRMM Operator multiplication with triangular matrix
%  FST_DTRMM(B,A,LEFT,TRS,{ALPHA=1})
%
%  B = ALPHA*op(A)*B (LEFT==1), op(X) == X or X', dep. on TRS
%  B = ALPHA*B*op(A) (LEFT==0)
% 
%  Here, A must be triangular (UPLO str. code required, DIAG code
%  optional) and square. If B has str. codes, they are ignored.
