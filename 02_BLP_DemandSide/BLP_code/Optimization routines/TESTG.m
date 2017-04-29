function f=testg(x)
% Usage:
% g=testg(x)
% TESTG is the interface function that used for
% accessing routines included to More set of test functions.
% TESTG returns the gradient <g> at a point <x>
%
global NDIM MDIM NPROB
f=grdfcn(NDIM,MDIM,x,NPROB); 
