function f=maxquadf(x)
% Usage:
% f=maxquadf(x)
% Returns the value <f> of Lemarechal's MAXQUAD function
% at a point <x>.
global matr_A matr_B mm
if size(x,2)>1, x=x;end
n=size(x,1);
f=(matr_A(1:n,:)*x)'*x - matr_B(:,1)'*x;
for l=2:mm
  d=(matr_A((l-1)*n+1:l*n,:)*x)'*x - matr_B(:,l)'*x;  if d>f, f=d; end
end
