%FST_ADDVEC Operator: add vector to cols/rows
%  FST_ADDVEC(C,X,COL,{ALPHA=1})
%
%  C = ALPHA*X*1' + C  (COL==1)
%  C = ALPHA*1*X' + C  (COL==0)
%
%  X vector, 1 vector of all ones. Structure codes ignored.
