function [xopt, fopt, nacc, nfcnev, nobds, ier, t, vm] = ...
    sim_anneal(fcn,x,max,rt,eps,ns,nt,neps,maxevl,lb,ub,c,iprint,t,vm,mymaxtime_sa);
%rand appears in lines 75, 84, 158
n = size(x,1);
xopt = zeros(n,1);
xp = zeros(n,1);
nacp = zeros(n,1);
tic_sa=cputime;
% *******************************************************************
% Set initial values
% *******************************************************************
nacc = 0;
nobds = 0;
nfcnev = 0;
ier = 99;
xopt = x;
nacp = zeros(n,1);
fstar = 1e20*ones(neps,1);
% *******************************************************************
% If the initial temperature is not positive, notify the user and abort
% *******************************************************************
 if (t <= 0.0);
   disp('The initial temperature is not positive. Reset the variable t');
   ier = 3; 
   return
end
% *******************************************************************
%  If the initial value is out of bounds, notify the user and abort
% *******************************************************************
 if (x>ub) & (x<lb)
   disp('initial condition out of bounds');
   ier = 2; r
   eturn
end
% *******************************************************************
% Evaluate the function with input x and return value as f.
% *******************************************************************
f = gmmobj(x);      
% *******************************************************************
%if the function is to be minimized, switch the sign of the function.
% Note that all intermediate and final output switches the sign back
%to eliminate any possible confusion for the user.
% *******************************************************************
if (max == 0) 
     f = -f;
end
% *******************************************************************
nfcnev = nfcnev + 1;
fopt = f;
fstar(1) = f;
% *******************************************************************
%   Start the main loop. Note that it terminates if (i) the algorithm
%   succesfully optimizes the function or (ii) there are too many
%   function evaluations (more than MAXEVL).
% *******************************************************************
quit=0;
while quit==0
    nup = 0;
    nrej = 0;
    nnew = 0;
    ndown = 0;
    lnobds = 0;
    
    m=1; 
    while m <= nt;
        j=1; 
        while j <= ns;
            h=1; 
            while h <= n;
% *******************************************************************            
% Generate xp, the trial value of x. Note use of vm to choose xp/
% *******************************************************************            
                i=1; 
                while i <= n;
                    if (i == h);
                        xp(i) = x(i) + (rand(1,1)*2.- 1.) * vm(i);
                    end
                    if (i ~= h);
                         xp(i) = x(i);
                    end
% *******************************************************************            
% If xp is out of bounds, select a point in bounds for the trial
% *******************************************************************            
                  if((xp(i) < lb(i)) | (xp(i) > ub(i)));
                      xp(i) = lb(i) + (ub(i) - lb(i))*rand(1,1);
                      lnobds = lnobds + 1;
                      nobds = nobds + 1;
                      if(iprint >= 3); 
                          disp('');
                          disp('current x');
                          disp(x)
                          if (max==1)
                             disp('current f');
                             disp(f);
                         end
                         if (max==0)
                             disp('current f');
                             disp(-f)
                         end
                         disp('trial x');
                         disp(xp)
                         disp('point rejected since out of bounds');
                      end
                  end
                  i = i+1; 
              end
% *******************************************************************            
%  Evaluate the function with the trial point xp and return as fp
% *******************************************************************
             fp = gmmobj(xp);    
             if(max == 0) 
                fp = -fp; 
             end
             nfcnev = nfcnev + 1;
% *******************************************************************            
% If too many function evaluations occur, terminate the algorithm.
% *******************************************************************            
             if(nfcnev >= maxevl);
                  disp('Too many function evaluations');
                  disp('consider increasing maxevl or eps');
                  disp('or decreasing nt or rt');
                  disp('These results are likely to be poor');
                      return
                  if (max == 0); 
                      fopt = -fopt; 
                  end
              ier = 1;
             end
% *******************************************************************                      
% Accept the new point if the function value increases
% *******************************************************************                      
              if (fp >= f);   %we have not closed it yet
                    if(iprint >= 3);
                        disp('point accepted');
                    end
                  x = xp;
                 f = fp;
                 nacc = nacc + 1;
                 nacp(h) = nacp(h) + 1;
                 nup = nup + 1;
             end
% *******************************************************************                      
% If greater than any other point, record as new optimum
% *******************************************************************                      
             if (fp > fopt)
               if(iprint >= 3); 
                   disp('new optimum'); 
               end
               xopt = xp;
               fopt = fp;
               nnew = nnew + 1;
             end
% *******************************************************************                      
%  if the point is lower, use the Metropolis criteria to decide on
%  acceptance or rejection.
% *******************************************************************                      
            if (fp<=fopt)
                 p = exp((fp - f)/t);
                 pp = rand(1,1);
                if (pp < p)
                    if(iprint >= 3); 
                        if (max)==1
                            disp('though lower, point accepted');
                        end
                        if(max)==0
                            disp('though higher, point accepted');
                        end                            
                    end
                     x = xp;
                     f = fp;
                     nacc = nacc + 1;
                     nacp(h) = nacp(h) + 1;
                     ndown = ndown + 1;
                end                     
                if (pp>=p)
                      nrej = nrej + 1;
                      if(iprint >= 3); 
                            if (max)==1
                                disp('lower point rejected');
                            end
                            if (max)==0
                               disp('higher point rejected');
                            end
                      end
                end
            end
            h = h+1; 
        end
    j = j+1; 
    end
% *******************************************************************                      
% Adjust vm so that approximately half of all evaluations are accepted.
% *******************************************************************                      
    i=1; 
    while i <= n;
        ratio = nacp(i) /ns;
        if (ratio > .6);
              vm(i) = vm(i)*(1. + c(i)*(ratio - .6)/.4);
          elseif (ratio < .4);
                 vm(i) = vm(i)/(1. + c(i)*((.4 - ratio)/.4));
        end
        if (vm(i) > (ub(i)-lb(i)));
               vm(i) = ub(i) - lb(i);
        end
        i = i+1; 
    end
 
   if(iprint >= 2);
      disp('-------------------------------------------------');       
      disp('intermediate results after step length adjustment');
      disp('-------------------------------------------------');
%       fprintf('elapsed cputime        : %18.6f\n',cputime-tic_sa);
%       if cputime-tic_sa>mymaxtime_sa
%           fprintf('Optimization terminated: time limit exceeded.\n');
%           return
%       end
      fprintf('new step vm            : \t');
      printm(vm');
      fprintf('current optimum theta2 : \t');
      printm(xopt');
      fprintf('current theta2         : \t');
      printm(x');
  end
   nacp = zeros(n,1);
   
m = m+1; 
end;
% *******************************************************************                      
if(iprint >= 1)
    totmov=nup+ndown+nrej;
        fprintf('intermediate results before next temperature reduction\n');
        fprintf('temperature:%8.3f\n',t);
    if max==1        
        fprintf('function value             :%18.4f\n',abs(fopt));
        fprintf('total moves                :%18.0f\n',totmov);
        fprintf('downhill                   :%18.0f\n',nup);
        fprintf('accepted uphill            :%18.0f\n',ndown);
        fprintf('rejected uphill            :%18.0f\n',nrej);
        fprintf('trials out of bounds       :%18.0f\n',lnobds);
        fprintf('new minima at this temp    :%18.0f\n',nnew);  
    end 
    if max==0
        fprintf('function value             :%18.4f\n',abs(fopt));
        fprintf('total moves                :%18.0f\n',totmov);
        fprintf('downhill                   :%18.0f\n',nup);
        fprintf('accepted uphill            :%18.0f\n',ndown);
        fprintf('rejected uphill            :%18.0f\n',nrej);
        fprintf('trials out of bounds       :%18.0f\n',lnobds);
        fprintf('new minima at this temp    :%18.0f\n',nnew);        
        
    end 
end
% *******************************************************************                      
% Check termination criteria
% *******************************************************************                       
 quit = 0;
 fstar(1) = f;
 if ((fopt - fstar(1)) <= eps) 
     quit = 1; 
 end
 if ( sum(abs(f-fstar)> eps) > 0 ) 
     quit = 0; 
 end;
% *******************************************************************                       
% Terminate SA if appropriate.
% *******************************************************************                       
if quit==1
   x = xopt;
   ier = 0;
   if (max == 0) 
      fopt = -fopt; 
   end
   if(iprint >= 1) 
     disp('SA achieved termination criteria. ier=0');
   end
end
% *******************************************************************                       
% If termination criteria are not met, prepare for another loop.
% *******************************************************************                       
t = rt*t;                            %temperature update
i=neps; 
while i >= 2;
  fstar(i) = fstar(i-1);
    i = i-1; 
end;                                            
f = fopt;                                       
x = xopt;                                      
end
