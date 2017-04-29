function f=illclinf(x)
% Usage:
% f=illclinf(x)
% Returns value <f> of the exact penalty function
% at a point <x> for the Ill-conditioned Linear Programing problem.
global matrA vectB vectC
if size(x,2)>1, x=x';end
n=size(x,1);
f=vectC'*x + 2*n*max([0;matrA*x-vectB;-x]);
