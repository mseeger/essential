%FST_DSPAMMT2 Operator multiplication by A*B', A, B sparse
%  FST_DSPAMMT2(D,A,B,C,{ALPHA=1},{BETA=0})
%
%  D = ALPHA*A*B'*C + BETA*D
% 
%  Here, A, B are sparse matrices with the same number of cols, while
%  D, C are dense matrices. Structure codes are ignored. Intermed.
%  storage of size COLS(A) is used, but B'*C does not have to be
%  stored.
