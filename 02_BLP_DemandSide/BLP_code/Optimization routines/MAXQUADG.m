function g=maxquadg(x)
% Usage:
% g=maxquadg(x)
% Returns the gradient <g> of Lemarechal's MAXQUAD function
% at a point <x>.
global matr_A matr_B mm
if size(x,2)>1, x=x;end
n=size(x,1);
f=(matr_A(1:n,:)*x)'*x - matr_B(:,1)'*x; k=1;
for l=2:mm
  d=(matr_A((l-1)*n+1:l*n,:)*x)'*x - matr_B(:,l)'*x;  if d>f, f=d; k=l; end
end
g=2*matr_A((k-1)*n+1:k*n,:)*x - matr_B(:,k);
