%FST_ACCUMULATE Special accumulation function
%  RET = FST_ACCUMULATE(A,IV,B,{BTR})
%
%  Specialized accumulation function. A is NN-by-NC matrix, IV an
%  NN index vector, B cell array of index vectors. For each i,
%  RET(i) is the sum of A(i,j), where j runs over B{IV(i)}.
%  If BTR is given, it is an additional index. In this case, RET(i)
%  is the sum of A(i,BTR(j)), where j runs over B{IV(i)}.
