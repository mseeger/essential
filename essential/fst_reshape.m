%FST_RESHAPE Operator for RESHAPE
%  FST_RESHAPE(A,M,N)
%
%  Replaces A by RESHAPE(A,M,N). The buffer of A has to have size
%  >= M*N.
%  ATTENTION, UNSAFE COMMAND:
%  Does not work if A is a field in a structure. Problem is not
%  resolved!
