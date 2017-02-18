function f = lmda(beta)
% This function computes the weights for the expected maximum log likelihood
global n N nb1d nb1s nb wght lambda pn
load 'Data/jecnew'; 

%% Set the value of GM, SG, B
GM=[1 -beta(nb1d+1+nb1s+1);-beta(nb1d+1) 1];
SG=[(beta(nb-2)) beta(nb);beta(nb) (beta(nb-1))];
B=[beta(1:nb1d)' 0 0 0 0 0; beta(nb1d+1+1) 0 beta(nb1d+1+2:nb1d+1+nb1s)']';

%% Update series
W_1=W; W_1(:,size(W,2))=ones(T,1);  %po replaced with ones
W_0=W; W_0(:,size(W,2))=zeros(T,1); %po replaced with zeros
wght_1=wght; %previous iteration weight
lambda_1=lambda; %previous period lambda
h1=diag(exp(-.5.*((Y*GM-W_1*B)*inv(SG)*(Y*GM-W_1*B)')));
h0=diag(exp(-.5.*((Y*GM-W_0*B)*inv(SG)*(Y*GM-W_0*B)')));
wght=(lambda_1*h1)./(lambda_1*h1+(1-lambda_1)*h0);
f=mean(wght);
tmp=corrcoef(wght,wght_1);
if tmp(1,2)>.999;
    n=N+1;
    pn=wght>.5;
else
    pn = wght;
end
