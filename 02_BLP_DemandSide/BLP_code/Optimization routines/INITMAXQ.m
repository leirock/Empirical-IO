function x=initmaxq(n,m)
% Usage:
% x=initmaxq(n,m)
% Sets the global constant matrices dependent on the dimension (n,m)
% and returns the standard starting point <x>
% for the Lemarechal's MAXQUAD function.
global matr_A matr_B mm
mm=m;
for l=1:m  
  for i=1:n,
    for j=i+1:n
       matr_A(i+(l-1)*n,j)=exp(i/j)*cos(i*j)*sin(l); 
       matr_A(j+(l-1)*n,i)=matr_A(i+(l-1)*n,j);
    end
    matr_B(i,l)=exp(i/l)*sin(i*l);
  end
  for i=1:n, j=find(1:n~=i);
    matr_A(i+(l-1)*n,i)=i*abs(sin(l))/n + sum(abs(matr_A(i+(l-1)*n,j)));
  end
end
  for i=1:n, x(i)=1; end, x=x';
