%FST_FLIPDIMS Operator for flipping dimensions
%  FST_FLIPDIMS(A,N)
%
%  Operates on A, which is of size N*D1-by-D2. The flat buffer of
%  A is taken as N-by-D1-by-D2 tensor, and dimensions 2 and 3
%  are exchanged. The resulting matrix is reshaped to have D1
%  columns. In-place variant of
%    RESHAPE(PERMUTE(RESHAPE(A,N,D1,D2),[1 3 2]),N*D2,D1),
%  but Matlab PERMUTE is not in-place.
%
%  NOTE: This is not straightforward if D1~=D2. We use CACM Algorithm
%  467 which does transposition in place:
%    Brenner: Matrix Transposition in Place, CACM 16(11), 1973
%    Available: portal.acm.org
%  See also matrix/TransInPlace.h of STATSIM.
