%FST_DSPAMM Operator multiplication by sparse matrix
%  FST_DSPAMM(C,A,B,TRA,{ALPHA=1},{BETA=0})
%
%  C = ALPHA*op(A)*B + BETA*C, op(A)==A / A', dep. on TRA
% 
%  Here, A is a sparse matrix, while C, B are dense matrices.
%  Structure codes are ignored.
