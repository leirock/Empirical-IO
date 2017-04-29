**************************************************************************
*To accompany Knittel and Metaxoglou (2008)
* Estimation of Random Coefficient Demand Models: 
* Challenges, Difficulties and Warnings
* Knittel      : crknittel@ucdavis.edu
* Metaxoglou   : konstantinos.metaxoglou@bateswhite.com
***************************************************************************
clear 
set memo 2000m
set more off
capture log close

global csv_file "E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\BLP_zero_code\Optimization results\Optimization results.csv"
global out_file "E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\BLP_zero_code\Optimization results\Optimization results.dta"
#delimit ;

insheet 

optmethod
stvalue       
fcnevals       
exitinfo       
toc
fval           

price_mean   
const_mean
hpwt_mean
air_mean
mpd_mean
space_mean

price_sigma
const_sigma
hpwt_sigma
air_sigma
mpd_sigma

price_mean_se   
const_mean_se
hpwt_mean_se
air_mean_se
mpd_mean_se
space_mean_se

price_sigma_se
const_sigma_se
hpwt_sigma_se
air_sigma_se
mpd_sigma_se


gradients1_price_sigma
gradients1_const_sigma
gradients1_hpwt_sigma
gradients1_air_sigma
gradients1_mpd_sigma
gradients1_norm_inf

gradients2_price_sigma
gradients2_const_sigma
gradients2_hpwt_sigma
gradients2_air_sigma
gradients2_mpd_sigma
gradients2_norm_inf

gradients3_price_sigma
gradients3_const_sigma
gradients3_hpwt_sigma
gradients3_air_sigma
gradients3_mpd_sigma
gradients3_norm_inf

hessians_eig1  
hessians_eig2  
hessians_eig3  
hessians_eig4  
hessians_eig5  

hessians2_eig1 
hessians2_eig2 
hessians2_eig3 
hessians2_eig4 
hessians2_eig5 

using "E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\BLP_zero_code\Optimization results\Optimization results.csv", clear;

local myvarlist
price_mean   
const_mean
hpwt_mean
air_mean
mpd_mean
space_mean

price_sigma
const_sigma
hpwt_sigma
air_sigma
mpd_sigma

price_mean_se   
const_mean_se
hpwt_mean_se
air_mean_se
mpd_mean_se
space_mean_se

price_sigma_se
const_sigma_se
hpwt_sigma_se
air_sigma_se
mpd_sigma_se


gradients1_price_sigma
gradients1_const_sigma
gradients1_hpwt_sigma
gradients1_air_sigma
gradients1_mpd_sigma
gradients1_norm_inf

gradients2_price_sigma
gradients2_const_sigma
gradients2_hpwt_sigma
gradients2_air_sigma
gradients2_mpd_sigma
gradients2_norm_inf

gradients3_price_sigma
gradients3_const_sigma
gradients3_hpwt_sigma
gradients3_air_sigma
gradients3_mpd_sigma
gradients3_norm_inf

hessians_eig1  
hessians_eig2  
hessians_eig3  
hessians_eig4  
hessians_eig5  

hessians2_eig1 
hessians2_eig2 
hessians2_eig3 
hessians2_eig4 
hessians2_eig5;


#delimit cr

*deal with irregular values
foreach v of local myvarlist {
	qui capture replace `v'="" if `v'=="NaN"	
	*imaginary numbers
	qui capture replace `v'="" if index(`v',"i")
	qui capture destring `v', replace
}

foreach v of local myvarlist {
	qui capture replace `v'=. if fval==10000000000 
}

sort optmethod stvalue
qui gen optmethod_str=""
qui replace optmethod_str="QNewton 1"        if optmethod==1
qui replace optmethod_str="Simplex"          if optmethod==2
qui replace optmethod_str="Solvopt"          if optmethod==3
qui replace optmethod_str="Conjgrad"         if optmethod==4
qui replace optmethod_str="QNewton 2"        if optmethod==5
qui replace optmethod_str="GA JBES"          if optmethod==6
qui replace optmethod_str="SA"  	     if optmethod==7
qui replace optmethod_str="MADS"             if optmethod==8
qui replace optmethod_str="GPS"              if optmethod==9
qui replace optmethod_str="GA Matlab"        if optmethod==10

format *mean* *se* *grad* *hess* %15.4fc

order optmethod_str optmethod stvalue fcnevals exitinfo toc fval

sort optmethod stvalue

compress

save "E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\BLP_zero_code\Optimization results\Optimization results.dta", replace

*eof
