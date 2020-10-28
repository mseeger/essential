%FST_INVCHOL Compute inverse from Cholesky factor (in place)
%  INFO = FST_INVCHOL(FACT)
%
%  Overwrites Cholesky factor FACT of matrix A by inverse of A.
%  The factor is given in FACT, which is upper or lower triangular
%  (UPLO field required). The inverse overwrites FACT, but just the
%  triangle occupied by the factor is accessed. INFO is the return
%  code by LAPACK DPOTRI (0: success).
