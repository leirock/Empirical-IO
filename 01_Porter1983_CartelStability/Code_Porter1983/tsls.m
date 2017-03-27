% This file mimicks column one in Table 3 of the Porter (1983) paper
clear
format short g
format compact
load('Data/jec.txt');

%% Variables
T = size(jec,1); %Total number of weeks
t = jec(:,1);L = jec(:,3);Q = jec(:,4);po = jec(:,5);gr = jec(:,6);
lnQ = log(Q); %log sales
lngr = log(gr); %log prices
const = ones(T,1); %A vector of ones for the constant

% Months dummies
m1 =(t>=1&t<=4)|(t>=53&t<=56)|(t>=105&t<=108)|(t>=157&t<=160)|(t>=209&t<=212)|(t>=261&t<=264)|(t>=313&t<=316);
m2 =(t>=5&t<=8)|(t>=57&t<=60)|(t>=109&t<=112)|(t>=161&t<=164)|(t>=213&t<=216)|(t>=265&t<=268)|(t>=317&t<=320);
m3 =(t>=9&t<=12)|(t>=61&t<=64)|(t>=113&t<=116)|(t>=165&t<=168)|(t>=217&t<=220)|(t>=269&t<=272)|(t>=321&t<=324);
m4 =(t>=13&t<=16)|(t>=65&t<=68)|(t>=117&t<=120)|(t>=169&t<=172)|(t>=221&t<=224)|(t>=273&t<=276)|(t>=325&t<=328);
m5 =(t>=17&t<=20)|(t>=69&t<=72)|(t>=121&t<=124)|(t>=173&t<=176)|(t>=225&t<=228)|(t>=277&t<=280);
m6 =(t>=21&t<=24)|(t>=73&t<=76)|(t>=125&t<=128)|(t>=177&t<=180)|(t>=229&t<=232)|(t>=281&t<=284);
m7 =(t>=25&t<=28)|(t>=77&t<=80)|(t>=129&t<=132)|(t>=181&t<=184)|(t>=233&t<=236)|(t>=285&t<=288);
m8 =(t>=29&t<=32)|(t>=81&t<=84)|(t>=133&t<=136)|(t>=185&t<=188)|(t>=237&t<=240)|(t>=289&t<=292);
m9 =(t>=33&t<=36)|(t>=85&t<=88)|(t>=137&t<=140)|(t>=189&t<=192)|(t>=241&t<=244)|(t>=293&t<=296);
m10=(t>=37&t<=40)|(t>=89&t<=92)|(t>=141&t<=144)|(t>=193&t<=196)|(t>=245&t<=248)|(t>=297&t<=300);
m11=(t>=41&t<=44)|(t>=93&t<=96)|(t>=145&t<=148)|(t>=197&t<=200)|(t>=249&t<=252)|(t>=301&t<=304);
m12=(t>=45&t<=48)|(t>=97&t<=100)|(t>=149&t<=152)|(t>=201&t<=204)|(t>=253&t<=256)|(t>=305&t<=308);
month=[m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12];

% Structural dummies
DM1 = t>=28 & t<=166;
DM2 = t>=167 & t<=181;
DM3 = t>=182 & t<=323;
DM4 = t>=324 & t<=328;
DM=[DM1,DM2,DM3,DM4];

% Regressors and instruments
X1=[const,L,month]; %Demand exogenous variables
X2=[const,DM,po,month]; %Supply exogenous variables
Z1=[X1,DM,po];
Z2=[X2,L];
X1=[X1,lngr]; 
X2=[X2,lnQ];

%% 2SLS

% Demand
W1=inv(Z1'*Z1); 
X1hat=Z1*W1*Z1'*X1;
beta1=inv(X1hat'*X1hat)*X1hat'*lnQ; %Estimated demand parameters

uhat1=lnQ-X1*beta1; %Demand residuals
shat1=(sum(uhat1.^2))/(T-size(X1,2)); %Estimated varience
se1=sqrt(diag(shat1*(inv(X1hat'*X1hat)))); %Estimated standard error
lnQhat=X1*beta1;

% Supply
W2=inv(Z2'*Z2);
X2hat=Z2*W2*Z2'*X2;
beta2=inv(X2hat'*X2hat)*X2hat'*lngr; %Estimated supply parameters

uhat2=lngr-X2*beta2; %Supply residuals
shat2=(sum(uhat2.^2))/(T-size(X2,2)); %Estimated varience
se2=sqrt(diag(shat2*(inv(X2hat'*X2hat)))); %Estimated standard error
lngrhat=X2*beta2;

shat12=(sum(uhat2.*uhat1))/(T-size(X2,2)); %Estimated covarience

param=[beta1;beta2(1);beta2(7:18);beta2(2:6);beta2(19);shat1;shat2;shat12];
save 'Data/param' param beta1 beta2 %They will be used as starting values later in MLE

R2_1=((sum((lnQ-mean(lnQ)).*(lnQhat-mean(lnQhat))))^2)/(sum((lnQ-mean(lnQ)).^2)*sum((lnQhat-mean(lnQhat)).^2));
R2_2=((sum((lngr-mean(lngr)).*(lngrhat-mean(lngrhat))))^2)/(sum((lngr-mean(lngr)).^2)*sum((lngrhat-mean(lngrhat)).^2));

%% Report estimated results
horz='Demand';
vert=['C      ';
      'LAKES  ';
      'GR     ';];
disp(['       ' horz])
for i=1:size(vert,1)
    if i<=size(vert,1)-1
        disp([vert(i,:),num2str(beta1(i))])
        disp(['       ','(',num2str(se1(i)),')'])
    else
        disp([vert(i,:),num2str(beta1(size(X1,2)))])
        disp(['       ','(',num2str(se1(size(X1,2))),')'])
    end
end
disp(['R2_1  ',num2str(R2_1)])
disp(['std1  ',num2str(sqrt(shat1))])

disp('       ')

horz='Supply';
vert=['C      ';
      'DM1    ';
      'DM2    ';
      'DM3    ';
      'DM4    ';
      'PO     ';
      'TQG    ';];
disp(['       ' horz])
for i=1:size(vert,1)
    if i<=size(vert,1)-1
        disp([vert(i,:),num2str(beta2(i))])
        disp(['       ','(',num2str(se2(i)),')'])
    else
        disp([vert(i,:),num2str(beta2(size(X2,2)))])
        disp(['       ','(',num2str(se2(size(X2,2))),')'])
    end
end
disp(['R2_2   ',num2str(R2_2)])
disp(['std2   ',num2str(sqrt(shat2))])