function g = apprgrdn(x,f,fun,deltax,obj)
% Usage:
% g = apprgrdn(x,f,fun,deltax,obj)
% Function apprgrdn.m performs the finite difference approximation 
% of the gradient <g> at a point <x>.
% <f> is the calculated function value at a point <x>,
% <fun> is the name of the Matlab function, which calculates function values
% <deltax> is a vector of the relative stepsizes,
% <obj> is the flag indicating whether the gradient of the objective
%        function (1) or the constraint function (0) is to be calculated. 
%
      n=max(size(x)); ee=ones(size(x));
        di=abs(x); idx=find(di<5e-15); di(idx)=5e-15*ee(idx); 
        di=deltax.*di; 
        if obj, idx=find(abs(di)<2e-10); di(idx)=2e-10*sign(di(idx));
        else,   idx=find(abs(di)<5e-15); di(idx)=5e-15*sign(di(idx));
        end
        y=x; 
        for i=1:n
          y(i)=x(i)+di(i);
          fi=feval(fun,y);
          if obj, if fi==f, 
             for j=1:3
                di(i)=di(i)*10;  y(i)=x(i)+di(i); 
                fi=feval(fun,y); if fi~=f, break, end
             end   
          end, end
          g(i)=(fi-f)/di(i);
          if obj, if any(idx==i)
            y(i)=x(i)-di(i);
            fi=feval(fun,y);
            g(i)=.5*(g(i)+(f-fi)/di(i));
          end, end            
          y(i)=x(i);
        end
