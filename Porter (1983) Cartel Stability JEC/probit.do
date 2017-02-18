clear
cd /Users/apple/Dropbox/Documents/Resources/06_Industrial_Organization/04_CU_XiaoEmprical_IO/Assignment1/Codes

use Data/jecnew
gen trend=_n
probit po L DM1-DM4 m1-m11
predict phat
replace phat=.5 if phat==. /* 5 cases, to avoid missing values in matlab, it will not affect the results */
sort trend
keep phat

save Data/wght01.csv, nol  replace
