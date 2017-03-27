% function [iter,bf0,bestfun,inarow,psi,beta,stopcode]=gamin(funstr,parspace,popseed,options)
function [ga_out,ga_beta,gafuneval]=gamin(funstr,parspace,popseed,options)

global mymaxtime mymaxfunevals

options.GenerationSize=100;

options.CrossoverRate=-1;     %If -1, use Booker's VCO, otherwise must be in (0,1).
options.MutationRate=0.02;    %If -1, use adaptive mutation rate, otherwise must be in (0,1).
options.DisplayIterations=5;  %Display status once every DisplayIterations generations.
                              %Set to zero to suppress printing.
options.MaximumIterations=20000;
options.Epsilon=1e-4;         %Smallest gain worth recognizing.
options.StopIterations=2000;  %Minimum number of gains < Epsilon before stop.
options.ReplaceBest=50;       %Every ReplaceBest generations, re-insert best-so-far 
options.VectorizedFunction=0; %Logical indicator set to 1 if the function can simultaneously 
options.ReproductionSelection='tournament';
options.RandomState=10;
% [ roulette | {tournament} | elitist ];
%   RandomState [seed based on clock time]

gafuneval=0;
m=options.GenerationSize;
eta=options.CrossoverRate;
gam=options.MutationRate;
printcnt=options.DisplayIterations;
printcnt=1;
maxiter=options.MaximumIterations;
epsln=options.Epsilon;
stopiter=options.StopIterations; 
rplcbest=options.ReplaceBest;
vecfun=options.VectorizedFunction;
repromethod=options.ReproductionSelection;
randomstate=options.RandomState;

% initialize random seed
rand('state',randomstate)

% Use Booker's VCO if eta==-1
vco=(eta<0);  eta=abs(eta);

% Use adaptive mutation rate if gam<0
adaptgamma=(gam<0);

% Cancel rplcbest if <=0 or if repromethod==elitist
if rplcbest<=0 | repromethod==2
  rplcbest=maxiter+1; 
end
K=size(parspace,1);

% Draw initial Generation
G=rand(K,m).*repmat(diff(parspace,1,2),1,m)+repmat(parspace(:,1),1,m);
G(:,1:size(popseed,2))=popseed;

% Initial 'best' holders
inarow=0;
bestfun=Inf; beta=zeros(K,1);
iter=0;  stopcode=0;
oldpsi=1;  % for VCO option

while stopcode==0
   iter=iter+1;
   
   % Locate unique chromosomes to minimize evaluations of funstr
   [Gcode,ix]=sort(rand(1,K)*G);   % 1xK vector
   G=G(:,ix);
   isfirst=find(([1 diff(Gcode)])>0);   % positions of unique chromosomes
   Gunq=G(:,isfirst);  munq=length(isfirst);
   
clonepos=sum(repmat((1:m)',1,length(isfirst))>=repmat(isfirst,m,1),2);
   
   % psi is number of unique rows in G divided by m.
   psi=munq/m;

   % Call function for each vector in G
   switch vecfun
     case 'on'
    for jj=1:size(Gunq,2);     
	    f(jj)=feval(funstr,Gunq(:,jj)); 
        best_fun00=bestfun;
        gafuneval=gafuneval+1
        if gafuneval==mymaxfunevals              
            best_fun00=bestfun;
            stopcode=1;
            break
        end
    end    
     otherwise
        f=zeros(1,munq);
        for i=1:munq
       	  f(i)=feval(funstr,Gunq(:,i));  
          best_fun00=bestfun;
         gafuneval=gafuneval+1;
         if gafuneval==mymaxfunevals
            best_fun00=bestfun;
            stopcode=1;
            break
        end

        end
   end

   % set NaN or infinite values to max(f)+1
   ee=(isnan(f) | isinf(f));
   if any(ee)
     f(ee)=repmat(max(f(~ee))+1,1,sum(ee));
   end
   clear ee

   % Restore full population of fitness scores
   f=f(clonepos);
   
   % Check for gains
   f0=f;
   [bf0,bx]=min(f);
   bf=min([bf0 bestfun]);
   fgain=(bestfun-bf);
   if fgain>epsln
	inarow=0;
   else
        inarow=inarow+1;
   end
   if fgain>0
	bestfun=bf;
    ga_beta=beta;
	beta=G(:,bx(1));
   end
   if printcnt>0 & (iter==0 | rem(iter,printcnt)==0)        
      fprintf('iteration  :%5i\t',iter);
%     fprintf('bf0        :%12.4f\n',bf0);
      fprintf('bestfun    :%12.4f\t',best_fun00);
      fprintf('inarow     :%5i\t',inarow);
      fprintf('100*psi    :%8.3f\t',100*psi);  
      fprintf('gafuneval  :%5i\n',gafuneval);   
%       ga_beta=beta;
   end
  if gafuneval==mymaxfunevals
      fprintf('GAMIN: Number of Function evaluations exceeded.\n');
      ga_out=[iter,abs(bf0),abs(best_fun00),inarow,100*psi];
%       ga_beta=beta;
      stopcode=1;
      break
  end

   % Reproduction
   switch repromethod
   case 'roulette'
     f1=(max(f)-f).^(1+log(iter)/100);
     if all(f1==0)
       pcum=((1:m)/m)';
     else
       pcum=cumsum(f1')/sum(f1);
     end
     r=rand(1,m); r=sum(repmat(r,m,1)>repmat(pcum,1,m))+1;
     clear f1 pcum
   case 'tournament'
     r=ceil(rand(2,m)*m);
     [f1,rx]=min(f(r));
     r=r((0:m-1)*2+rx);
     clear f1 rx
   case 'elitist'
     [f1,rx]=sort(f);
     r=[rx(randperm(m/2)) rx(randperm(m/2))];
     clear f1 rx
   otherwise
     error('Invalid ReproductionSelection option.')
   end
   G=G(:,r);
   
   % Crossover
   if vco
        % psi is number of unique rows in G divided by m.
	eta=max(min(eta-psi+oldpsi,1),0.2);
	oldpsi=psi;
   end   
   y=sum(rand(m/2,1)<eta); 
   if y>0
   
       % choose crossover point
     x=ceil(rand(y,1)*(K-1));
     for i=1:y
        G(x(i)+1:K,i+[0 m/2])=G(x(i)+1:K,i+[m/2 0]);
     end
   end

   % Mutation
   if adaptgamma
     gam=0.005 + 0.1*(1-sqrt(psi));
   end
   
M=rand(K,m).*repmat(diff(parspace,1,2),1,m)+repmat(parspace(:,1),1,m);
   domuta=find(rand(K,m)<gam);
   G(domuta)=M(domuta);
   
   % Once every rplcbest generations, re-insert best beta
   if rem(iter,rplcbest)==0
        G(:,m)=beta;
   end
   stopcode=(inarow>stopiter)+2*(iter>maxiter);
end
stopcode=logical(stopcode-1);  % 0 if normal stop, 1 if maxiter exceeded

% if printcnt>0
%    if stopcode
%      disp(sprintf('GAMIN: Maximum number of generations exceeded.\n'))
%    else
%      disp(sprintf('GAMIN: No improvement in %5.0f generations.\n',stopiter))
%    end
% end
% ga_out=[iter -abs(bf0) -abs(best_fun00) inarow 100*psi];
return 
%% end of gamin.m

%[beta,stopcode]=GAMIN(funstr,parspace,popseed,options,p1,p2,p3,p4,...)
%
% Genetic Algorithm for function minimization.
% Especially useful for functions with kinks and discontinuities
% and where a good "starting point" is unavailable.
%    See Dorsey and Mayer, 
%    Journal of Business and Economic Statistics, January 1995, 13(1)
% Program by Michael Gordy <mgordy@frb.gov>
% Version 3.1, 26 January 2001, written for Matlab v5
%
% INPUTS:
%  funstr     = (String) Name of objective function to be minimized.
%               Function should map (Kx1) parameter vector to a scalar
%               output or, if vectorized, (KxL) parameter matrix to
%               a (1xL) vector of function values for the L parameter 
%               column vectors.
%  parspace   = (K x 2) matrix is [min max] of parameter space dimensions.
%  popseed    = (K x m1) matrix of starting values with which to seed
%               the initial population.  Need m1<=GenerationSize.
%               (optional).
%  options    = structure variable of option settings (optional).
%  p1,p2,...  = Additional data to be passed to funstr (optional).
%
% The fields of options are [default value in brackets]:
%   GenerationSize [20]
%           Must be even integer.
%   CrossoverRate [-1]
%           If -1, use Booker's VCO, otherwise must be in (0,1).
%   MutationRate [0.02]
%           If -1, use adaptive mutation rate, otherwise must be in (0,1).
%   DisplayIterations [0]
%           Display status once every DisplayIterations generations.
%           Set to zero to suppress printing.
%   MaximumIterations [20000]
%   Epsilon [1e-4]
%           Smallest gain worth recognizing.
%   StopIterations [2000]
%           Minimum number of gains < Epsilon before stop.
%   ReplaceBest [50]
%           Every ReplaceBest generations, re-insert best-so-far 
%           into population.
%   VectorizedFunction [0]
%           Logical indicator set to 1 if the function can simultaneously 
%           evaluate many parameter vectors.
%   ReproductionSelection [ roulette | {tournament} | elitist ];
%   RandomState [seed based on clock time]
%
% OUTPUTS:
%   beta       = (Kx1) Parameter vector minimizing funstr
%   stopcode   = Code for terminating condition
%                 == 0 if terminated normally
%                 == 1 if maximum number of generations exceeded
%
% NOTES: 
% (1) Objective function is minimized with respect to its first parameter,
%     which is expressed as a column vector.
%     Example: 
%       Say we want to minimize function f with respect to vector p,
%       and need also to pass to f data matrices x,y,z.  Then,
%       write the function f so it is called as f(p,x,y,z).  GA will
%       assume that p is a column vector.
% (2) Big gains in speed if objective function can be vectorized.
% (3) Automatic elimination of vectors which return function value of NaN.
% (4) Intermediate results are saved to "gabest.mat".  This allows
%     you to pick up where you left off after an interruption.