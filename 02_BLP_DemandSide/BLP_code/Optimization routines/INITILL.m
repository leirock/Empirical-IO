function x=initill(n)
% Usage:
% x=initill(n)
% Sets the global constant matrices and vectors of the dimension <n>
% and returns the standard starting point <x>
% for the Ill-conditioned Linear Programming problem.
global matrA vectB vectC
vectB=[]; vectC=[]; matrA=[];  
  for i=1:n,
    x(i)=0;     vectB(i)=0;
    for j=1:n
       matrA(i,j)=1/(i+j);       vectB(i)=vectB(i)+1/(i+j);
    end
    vectC(i)=-1/(i+1)-vectB(i);
  end
  vectB=vectB';    vectC=vectC';  x=x';
