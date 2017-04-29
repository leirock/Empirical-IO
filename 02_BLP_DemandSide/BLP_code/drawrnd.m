clear
clc

income=normrnd(18147,13936.32871,6000,50);
income=income/10000;
income2=income.^2;
hhsize=normrnd(2.3,0.232307,6000,50);
child=normrnd(1.75,0.20984121,6000,50);

demogr=[income,income2,hhsize,child];
v=normrnd(0,1,90,250);
save draw v demogr