function f=testf(x)
% Usage:
% f=testf(x)
% TESTF is the interface function that used for
% accessing routines included to More set of test functions.
% TESTF returns the function value <f> at a point <x>
%
global NDIM MDIM NPROB
f=objfcn(NDIM,MDIM,x,NPROB); 
