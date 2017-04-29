function f=shorf(x)
% Usage:
% f=shorf(x)   
% Returns the value <f> of Shor's piece-wise quadratic function
% at a point <x>.
a=[ 0,  2,  1,  1,  3,  0,  1,  1,  0,  1;...
    0,  1,  2,  4,  2,  2,  1,  0,  0,  1;...
    0,  1,  1,  1,  1,  1,  1,  1,  2,  2;...
    0,  1,  1,  2,  0,  0,  1,  2,  1,  0;...
    0,  3,  2,  2,  1,  1,  1,  1,  0,  0];
b=[1;5;10;2;4;3;1.7;2.5;6;4.5];
if size(x,2)>1
  x=x';
end
  f=0;
  for i=1:10;
     d=b(i)*sum((x-a(1:5,i)).^2);
     if d>f,f=d;end
  end
