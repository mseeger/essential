%FST_DSYRK Operator outer product
%  FST_DSYRK(C,A,TRS,{ALPHA=1},{BETA=0})
%
%  C = ALPHA*op(A)*op(A)' + BETA*C, op(X) == X or X', dep. on TRS
% 
%  Here, C must be symmetric (UPLO str. code required). A str. codes
%  are ignored.
