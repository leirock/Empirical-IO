clear
*cd /Users/apple/Dropbox/Documents/Academic/06_IO/04_CU_Emprical/hw2_BLP/Data
cd E:\Dropbox\Documents\Academic\06_IO\04_CU_Emprical\hw2_BLP\Data
use car_bel_UK

sort ye ma frm
egen cdid=group(ye ma)  //cdid
rename frm firm         //firm id
rename co mdl           //model id
rename princ price      //price
gen hpwt=hp/we          //hpwt
gen const=1             //constant term
gen space=le*wi*he

*qui tabulate brd, gen(brd) //brand dummy
gen share=qu/pop        //market share

*duplicates drop ye ma co,force
*gen co_str=string(co)
*replace co_str = "00" + co_str if length(co_str) == 1
*replace co_str = "0" + co_str if length(co_str) == 2
*gen id = real(string(ye) + string(ma) + co_str)
//gen id variable

keep cdid const firm hpwt home price share mdl space cy
*drop brd
aorder
*order brd*,last
export delimited using autoblp.csv, novarnames nolabel replace

**********************
** Generate cdindex **
**********************
bysort cdid: egen cdindex0=count(cdid)  //cdindex
duplicates drop cdid,force //cdindex for output

**********************
** Compare the blow **
**********************
*bysort cdid: egen aaa=sum(const)
*bysort cdid: egen bbb=count(const)
*bysort cdid: egen cdindex0=count(cdid)
*egen a=sum(cdindex0)
*gen a1=sum(cdindex0)
* collapse

************************
** Summary Statistics **
************************
estpost tabstat cy home hpwt  price share space, /* 
 */  statistics( mean sd min max ) /* 
 */  columns(variables)
esttab using myfile.tex,  /* 
 */ cells("cy home hpwt  price share space") /* 
 */  replace nonum noobs
