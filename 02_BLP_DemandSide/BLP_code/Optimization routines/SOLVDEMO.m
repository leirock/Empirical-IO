format long;
echo on; clc
%Tutorial For The Solver For Local Nonlinear Optimization Problems (SolvOpt)
%###########################################################################
%
%   This is a demonstration script-file for the SolvOpt
%   It closely follows the Tutorial section of the users' guide
%
%   It consists of:   
%     (1) minimization of a sample nonsmooth function (Problem 1),
%     (2) minimization of the user pre-defined penalty function
%         for constrained sample problem 2,
%     (3) solution of constrained minimization sample problem 3,
%     (4) illustration to the use of optional parameters to the solver.
%
%############################################################################
%
pause % Strike any key to continue (Use Ctrl-C to abort)
clc
%############################################################################
%       DEMO 1: UNCONSTRAINED PROBLEM
%---------------------------------------------------------------------------
%
% Consider the problem of finding the minimum to 
% Shor's piece-wise quadratic function:
% f(x)=max{f (x): i=1,...,m}
%           i
%
pause % Strike any key to continue

% Step 1. Code M-function that returns the objective function value at a point
type shorf;
pause % Strike any key to continue

%   Take a guess at the solution

x=[-1;1;-1;1;-1];

pause % Strike any key to start the optimization


[x,f,options]=solvopt(x,'shorf'); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of function evaluations was: 
options(10)
pause % Strike any key for the next demo
clc

% Minimize Shor's function 
% with user-supplied analytically calculated gradients

pause % Strike any key to continue

% Step 1. Code M-function that returns the objective function value at a point
% Step 1. Code M-function that returns the (sub)gradient at a point
type shorg;

pause % Strike any key to continue

%   Take a guess at the solution

x=[-1;1;-1;1;-1];

pause % Strike any key to start minimization


[x,f,options]=solvopt(x,'shorf','shorg'); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of function evaluations was: 
options(10)
%  The total number of gradient evaluations was: 
options(11)
pause % Strike any key for the next demo
clc

%###########################################################################
%       DEMO 2: MINIMIZATION OF A USER-DEFINED PENALTY FUNCTION
%---------------------------------------------------------------------------
%
%
% Consider  the Ill-conditioned Linear Programming Problem
% 
% minimize f(x)=C'x
%
% subject to
%         x>=0, Ax-B<=0
%
pause % Strike any key to continue
% The problem is equivalent to the minimization
% of the following exact penalty function
%
%         f(x)=C'x + R*max(0,Ax-B,-x),
%  where
%         R is the penalty coefficient.
%  The exact value of the penalty coefficient R=2n is known a priori.
%
pause % Strike any key to continue
clc
% Step 1. Code M-function that sets the values for the constant
%         matrices and vectors and returns the starting point:
pause % Strike any key to continue
type initill;
pause % Strike any key to continue
% Step 2. Code M-function that returns the penalty function value at a point
pause % Strike any key to continue
type illclinf;
pause % Strike any key to continue

%   Set the constant values and the starting point

x=initill(15);  

pause % Strike any key to start the optimization


[x,f,options]=solvopt(x,'illclinf'); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of function evaluations was: 
options(10)
pause % Strike any key for the next demo
clc

%###########################################################################
%       DEMO 3: CONSTRAINED MINIMIZATION SAMPLE PROBLEM 3
%---------------------------------------------------------------------------
%
% Consider Shell Dual Problem
%
%                      3
% minimize f(y,z)=2(d,y )+(Cy,y)-(b,z)
%
% subject to
%              y>=0, z>=0,
%                            2
%       g(y,z)=Az-2Cy-e-3(d,y )<=0
%
pause % Strike any key to continue
%

pause % Strike any key to continue
clc
% Step 1. Code M-function that sets the values for the constant
%         matrices and vectors and returns the starting point:
pause % Strike any key to continue
type initdual;
pause % Strike any key to continue
% Step 2. Code M-function that returns the objective function value at a point
pause % Strike any key to continue
type dsobjf;
pause % Strike any key to continue
% Step 3. Code M-function that returns the (sub)gradient of the objective
%         function at a point
pause % Strike any key to continue
type dsobjg;
pause % Strike any key to continue
% Step 4. Code M-function that returns the maximal residual for a set of
%         constraints at a point
pause % Strike any key to continue
type dscntf;
pause % Strike any key to continue
% Step 5. Code M-function that returns the (sub)gradient of the constraint
%         with the maximal residual at a point
pause % Strike any key to continue
type dscntg;
pause % Strike any key to continue

%   Set the constant values and the starting point

x=initdual;

pause % Strike any key to start the optimization


[x,f,options]=solvopt(x,'dsobjf','dsobjg',[],'dscntf','dscntg'); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of objective function evaluations was: 
options(10)
%  The total number of objective function gradient evaluations was: 
options(11)
%  The total number of constraints evaluations was: 
options(12)
%  The total number of constraint function gradient evaluations was: 
options(13)
pause % Strike any key for the next demo
clc

%###########################################################################
%       DEMO 4: OPTIONAL PARAMETER SETTINGS
%---------------------------------------------------------------------------
%       Solve Problem 1 at 
%       required relative error for the argument OPTIONS(2)=1e-6 
%                                     (the default value is 1e-4) and
%       required relative error for the function OPTIONS(3)=1e-8
%                                     (the default value is 1e-6).
%
pause % Strike any key to continue
clc

%   Set the optional parameters

options=soptions; options(2)=1e-6; options(3)=1e-8;

%   Initialize the problem data

x=[-1;1;-1;1;-1];

pause % Strike any key to start the optimization

[x,f,options]=solvopt(x,'shorf','shorg',options); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of function evaluations was: 
options(10)
%  The total number of gradient evaluations was: 
options(11)

pause % Strike any key to continue
clc
%       Solve Problem 3 at 
%       the upper bound admissible maximal residual OPTIONS(2)=1e-12
%                                        (the default value is 1e-8).
%
pause % Strike any key to continue

options=soptions; options(6)=1e-12;

pause % Strike any key to continue
%   Initialize the problem data

x=initdual;

[x,f,options]=solvopt(x,'dsobjf','dsobjg',options,'dscntf','dscntg'); 

%  The solver has found a solution at:
x
%  The function value at the solution is:
f
pause % Strike any key to continue
%  The total number of objective function evaluations was: 
options(10)
%  The total number of objective function gradient evaluations was: 
options(11)
%  The total number of constraints evaluations was: 
options(12)
%  The total number of constraint function gradient was: 
options(13)


%  Thank you! Bye!

echo off
