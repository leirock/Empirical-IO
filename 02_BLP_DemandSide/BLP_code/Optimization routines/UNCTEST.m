function out=unctest(appr,number,filename)
% Usage:
% unctest(appr,number,filename)
% <unctest.m> is used to run Mor`e standard tests on <solvopt.m>
% with the default optional parameters.
% The input argument <appr> controls the way the gradients are calculated:
%   if appr==0, analytically calculated gradients are used (default value),
%   if appr==1, gradients are calculated by finite differences.
% The input vector <number> specifies the range of numbers of the
% test functions. The valid numbers are in the range 1:26,28:34.
% At number==0 (default value), the online dialog is invoked.
% The input argument <filename> specifies the name of the two files with
% extentions <.txt> and <.tbl> to write the results to.  
% By default, the two files <unctest.tbl> and <unctest.txt> will be
% created in the current Matlab working directory.

global NDIM MDIM NPROB

function_strings=[...
'#  1.  ROSE      2    2   Rosenbrock                             ' ;...
'#  2.  FROTH     2    2   Freudenstein and Roth                  ' ;...
'#  3.  BADSCP    2    2   Powell Badly Scaled                    ' ;...
'#  4.  BADSCB    2    3   Brown Badly Scaled                     ' ;...
'#  5.  BEALE     2    3   Beale                                  ' ;...
'#  6.  JENSAM    2   10   Jennrich and Sampson                   ' ;...
'#  7.  HELIX     3    3   Helical Valley                         ' ;...
'#  8.  BARD      3   15   Bard                                   ' ;...
'#  9.  GAUSS     3   15   Gaussian                               ' ;...
'# 10.  MEYER     3   16   Meyer                                  ' ;...
'# 11.  GULF      3   10   Gulf Research and Development          ' ;...
'# 12.  BOX       3   10   Box 3-Dimensional                      ' ;...
'# 13.  SING      4    4   Powell Singular                        ' ;...
'# 14.  WOOD      4    6   Wood                                   ' ;...
'# 15.  KOWOSB    4   11   Kowalik and Osborne                    ' ;...
'# 16.  BD        4   20   Brown and Dennis                       ' ;...
'# 17.  OSB1      5   33   Osborne 1                              ' ;...
'# 18.  BIGGS     6   13   Biggs EXP6                             ' ;...
'# 19.  OSB2     11   65   Osborne 2                              ' ;...
'# 20.  WATSON  ( 9)  31   Watson                                 ' ;...
'# 21.  ROSEX   (10) (10)  Extended Rosenbrock                    ' ;...
'# 22.  SINGX   (10) (10)  Extended Powell Singular               ' ;...
'# 23.  PEN1    ( 4) ( 5)  Penalty I                              ' ;...
'# 24.  PEN2    ( 4) ( 8)  Penalty II                             ' ;...
'# 25.  VARDIM  (10) (12)  Variably Dimensioned                   ' ;...
'# 26.  TRIG    (10) (10)  Trigonometric                          ' ;...
'# 27.  ALMOST  (10) (10)  Brown Almost Linear                    ' ;... % not available
'# 28.  BV      (10) (10)  Discrete Boundary Value                ' ;...
'# 29.  IE      (10) (10)  Discrete Integral Equation             ' ;...
'# 30.  TRID    (10) (10)  Broyden Tridiagonal                    ' ;...
'# 31.  BAND    (10) (10)  Broyden Banded                         ' ;...
'# 32.  LIN     (10) (20)  Linear --- Full Rank                   ' ;...
'# 33.  LIN1    (10) (20)  Linear --- Rank 1                      ' ;...
'# 34.  LIN0    (10) (20)  Linear - Rank 1 with Zero Cols. & Rows ' ];

function_min_values=[...
 0           , ...%01
 48.98425    , ...%02    2 local minima
 0           , ...%03
 0           , ...%04
 0           , ...%05
 124.36218   , ...%06
 0           , ...%07
 8.21487e-03 , ...%08
 1.127932e-8 , ...%09
 87.9458     , ...%10
 0           , ...%11
 0           , ...%12    a number of local minima
 0           , ...%13
 0           , ...%14
 3.075055e-4 , ...%15    additionally takes local minima on a surface
 85822.2     , ...%16
 5.46489e-5  , ...%17    a number of areas of flatness, but the only one minimum
 0           , ...%18    at least 3 local minima
 4.01377e-2  , ...%19    a number of local minima /  f'(x(2:4,6:11))=0
 1.39976e-6  , ...%20
 0           , ...%21
 0           , ...%22
 2.249975e-5 , ...%23
 9.376293e-6 , ...%24
 0           , ...%25
 2.79506e-5  , ...%26     +2 local minima
 0           , ...%27
 0           , ...%28
 0           , ...%29
 0           , ...%30     a number of local minima (essentially multiextr.)
 0           , ...%31     at least 3 local minima
 10          , ...%32
 4.63415     , ...%33     additionally takes a minimum at a surface
 6.13514     ];   %34     additionally takes a minimum at a surface

format long;

if nargin==0,                appr=0;
elseif isempty(appr),        appr=0
elseif appr~=0 & appr~=1,    appr=0;
end
if nargin<2,    number=0;
elseif isempty(number) | all(number>34) | number==27
   disp('Numbers out of the range'); return
end

if nargin<3
   filename1='unctest.tbl';
   filename2='unctest.txt';
else,
   filename1=setstr([filename,'.tbl']);
   filename2=setstr([filename,'.txt']);
end
fd=fopen(filename1,'wt');
ft=fopen(filename2,'wt');

range=[1:26,28:34]; nt=34;

if number==0
 default=1;
 
  disp('###################  MOR`E SET OF TEST FUNCTIONS  ##################');
  disp(' ');

  % Dialog 1. Enter the number of a test function =========================

  k=0; number=-1; 
  while number==-1
  disp('Enter the number of a test function from the list below or 0 to stop');
    if k>=nt, k=0; end
    for  j=k+1:min(k+12,nt),  
         if j~=27, disp(function_strings(j,:));  end
    end     
    k=k+12;  disp('Press Enter for more...');
    s=input('>>> ','s'); 
    if s=='0',disp('Bye'); fclose(fd); fclose(ft); return, end
    if ~isempty(s)
       number=sscanf(s,'%i');
       if ~any(range==number),number=-1; disp('Not in the range!'); end
    end  
  end
  
  % Initialize a problem data ==============================================
  
  NPROB=number;
  [NDIM,MDIM,x0]=initf(NPROB);
  
  % Dialog 2. Enter a starting point =======================================

   disp('The standard starting point is'); 
   disp(sprintf(' %g;',x0)); 
   disp('To accept it press Enter');
   disp(sprintf...
   ('To start at another one enter the %i coordinates separated by blanks',...
     NDIM));
   while 1,
     s=input('>>> ','s');
     if isempty(s), break
     elseif size(s,2)==NDIM, x0=sscanf(s,'%g'); break
     else,  disp('Wrong number of coordinates. Repeat entering');
     end
   end

  % Dialog 3. Enter a starting point =======================================

   disp('Would you like to use analytically calculated gradients? [y]'); 
   while 1,
     s=input('>>> ','s');
     if isempty(s), appr=0, break
     elseif s=='y' | s=='Y', appr=0; break
     elseif s=='n' | s=='N', appr=1; break 
     else, disp('Please, enter "y" or "n"'); 
     end
   end

%====================================================================

else,  default=0;
 nt=max(size(number)); if nt>1 & nt==size(number,1), number=number'; end
 i=1; while i<=nt, 
        if ~any(range==number(i)),  
           number=[number(1:i-1),number(i+1:nt)]; nt=nt-1; 
        else,  i=i+1;
        end
      end
 if nt==0, 
  disp('Numbers out of the range'); fclose(fd); fclose(ft); return
 end
end

if  appr 
fprintf(fd,'\n      FMinValue      FTrueMinValue  |f-f*|/|f|  FunEvaluatn');
else
fprintf(fd,'\n      FMinValue      FTrueMinValue  |f-f*|/|f|  Fnctn Grdnt');
end


nnf=0; nng=0;

for NPROB=number
    
  if ~default, [NDIM,MDIM,x0]=initf(NPROB); end
  fun='testf';grad='testg';
  x=x0; options=soptions;

  disp(sprintf('\n\n') );          fprintf(ft,'\n\n');
  disp(function_strings(NPROB,:)); fprintf(ft,'\n%s',function_strings(NPROB,:));
  disp('Starting Point is:');      fprintf(ft,'\n Starting point:%a');
  disp(sprintf('%g; ',x0) );       fprintf(ft,'%g; ',x0);

  if appr,       [x,f,options]=solvopt(x,fun,[],options);
  else,          [x,f,options]=solvopt(x,fun,grad,options);
  end

  if options(9) < 0,  fprintf(ft,'\nABNORMAL TERMINATION. CODE = %i.',options(9));
  else,               fprintf(ft,'\nNORMAL TERMINATION');    end    
   
   fprintf('\nValue of the function at the solution: %22.15g', f);
   fprintf(ft,'\nValue of the function at the solution: %22.15g', f);
   fprintf('\nNumber of function evaluations: %i', options(10));
   fprintf(ft,'\nNumber of function evaluations: %i', options(10));
   if ~appr
   fprintf('\nNumber of gradient evaluations: %i', options(11));
   fprintf(ft,'\nNumber of gradient evaluations: %i', options(11));
   end
   fprintf('\nMinimum Point x:');
   fprintf(ft,'\nMinimum Point x:');
   fprintf('%22.15g;  ', x);
   fprintf(ft,'%22.15g;  ', x);

   fmin=function_min_values(NPROB);

%  Special cases
   if     NPROB==2,
       fmin02=0;         if abs(f-fmin)>abs(f-fmin02), fmin=fmin02;  end   
   elseif NPROB==6,
       fmin06=259.58019; if abs(f-fmin)>abs(f-fmin06), fmin=fmin06;  end   
   elseif NPROB==11,
       fmin11=[.038,.038,.0385,fmin];
       [df,index]=min(abs(fmin11-[f,f,f,f]));  
       fmin=fmin11(index);
   elseif NPROB==15,
       fmin15=[0.00179453906640,fmin];
       [df,index]=min(abs(fmin15-[f,f])); 
       fmin=fmin15(index);
   elseif NPROB==18,
       fmin18=[5.65565e-003,0.30636677262479,fmin]; 
       [df,index]=min(abs(fmin18-[f,f,f]));  
       fmin=fmin18(index);
   elseif NPROB==19,
       fmin19=[1.78981358688109,26.305657,fmin]; 
       [df,index]=min(abs(fmin19-[f,f,f]));  
       fmin=fmin19(index);
   elseif NPROB==26,
       fmin26=0;          if abs(f-fmin)>abs(f-fmin26), fmin=fmin26;  end   
   elseif NPROB==30,
       fmin30=[1.02865203567795,1.36025590473840,1.34953612374127,1.05122618838356,0.71260601731262,0.39737346895853,fmin];
       [df,index]=min(abs(fmin30-[f,f,f,f,f,f,f]));  
       fmin=fmin30(index);
   elseif NPROB==31,
       fmin31=[3.05727843,2.68021992072616,fmin];
       [df,index]=min(abs(fmin31-[f,f,f]));  
       fmin=fmin31(index);
   end

   if abs(fmin)<eps,  df=abs(f-fmin);
   else,              df=abs((f-fmin)/fmin);
   end
   
  
   nnf=nnf+options(10); nng=nng+options(11);
   
   if     ~appr
    fprintf(fd,...
    '\n%2i:  %13.5e  %13.5e  %13.5e  %5i %5i',...
    NPROB,f,fmin,df,options(10),options(11));
   
   else
    fprintf(fd,...
    '\n%2i:  %13.5e  %13.5e  %13.5e  %8i',...
    NPROB,f,fmin,df,options(10));
   
   end

end

fprintf(fd,'\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
if ~appr, fprintf(fd,'\n Total gradient evaluations = %i',nng); end
fprintf(fd,'\n Total function evaluations = %i',nnf);

fclose(fd); fclose(ft);
