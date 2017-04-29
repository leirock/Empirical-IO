function g=illcling(x)
% Usage:
% g=illcling(x) 
% Returns the gradient <g> of the exact penalty function
% at a point <x> for the Ill-conditioned Linear Programming problem.
global matrA vectB vectC
if size(x,2)>1, x=x';end
n=size(x,1);
g=vectC;
f=0;
for i=1:n
  d=matrA(i,:)*x-vectB(i); if d>f, f=d; k=i; end
end
for i=1:n
  if -x(i)>f, f=-x(i); k=i+n; end
end
if f>0,
  if k>n, g(k-n)=g(k-n)-2*n;
  else, g=g+2*n*matrA(:,k);
  end
end
