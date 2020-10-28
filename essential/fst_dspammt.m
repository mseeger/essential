%FST_DSPAMMT Operator multiplication by A*A', A sparse
%  FST_DSPAMMT(C,A,B,{ALPHA=1},{BETA=0})
%
%  C = ALPHA*A*A'*B + BETA*C
% 
%  Here, A is a sparse matrix, while C, B are dense matrices.
%  Structure codes are ignored. Intermed. storage of size COLS(B)
%  is used, but A'*B does not have to be stored.
