clear
cd /Users/apple/Dropbox/Documents/Resources/06_Industrial_Organization/04_CU_XiaoEmprical_IO/Assignment1/solution/Codes

import delimited Data/jec.txt, delimiter(space) 

// Rename variables
loc s=0
foreach i in t v2 L Q po gr{
loc s=`s'+1
ren v`s' `i'
}
drop v2

// Endogenous variables
gen lnQ = log(Q)
gen lngr = log(gr)

// Months dummies
gen m1 =(t>=1&t<=4)|(t>=53&t<=56)|(t>=105&t<=108)|(t>=157&t<=160)|(t>=209&t<=212)|(t>=261&t<=264)|(t>=313&t<=316)
gen m2 =(t>=5&t<=8)|(t>=57&t<=60)|(t>=109&t<=112)|(t>=161&t<=164)|(t>=213&t<=216)|(t>=265&t<=268)|(t>=317&t<=320)
gen m3 =(t>=9&t<=12)|(t>=61&t<=64)|(t>=113&t<=116)|(t>=165&t<=168)|(t>=217&t<=220)|(t>=269&t<=272)|(t>=321&t<=324)
gen m4 =(t>=13&t<=16)|(t>=65&t<=68)|(t>=117&t<=120)|(t>=169&t<=172)|(t>=221&t<=224)|(t>=273&t<=276)|(t>=325&t<=328)
gen m5 =(t>=17&t<=20)|(t>=69&t<=72)|(t>=121&t<=124)|(t>=173&t<=176)|(t>=225&t<=228)|(t>=277&t<=280)
gen m6 =(t>=21&t<=24)|(t>=73&t<=76)|(t>=125&t<=128)|(t>=177&t<=180)|(t>=229&t<=232)|(t>=281&t<=284)
gen m7 =(t>=25&t<=28)|(t>=77&t<=80)|(t>=129&t<=132)|(t>=181&t<=184)|(t>=233&t<=236)|(t>=285&t<=288)
gen m8 =(t>=29&t<=32)|(t>=81&t<=84)|(t>=133&t<=136)|(t>=185&t<=188)|(t>=237&t<=240)|(t>=289&t<=292)
gen m9 =(t>=33&t<=36)|(t>=85&t<=88)|(t>=137&t<=140)|(t>=189&t<=192)|(t>=241&t<=244)|(t>=293&t<=296)
gen m10=(t>=37&t<=40)|(t>=89&t<=92)|(t>=141&t<=144)|(t>=193&t<=196)|(t>=245&t<=248)|(t>=297&t<=300)
gen m11=(t>=41&t<=44)|(t>=93&t<=96)|(t>=145&t<=148)|(t>=197&t<=200)|(t>=249&t<=252)|(t>=301&t<=304)
gen m12=(t>=45&t<=48)|(t>=97&t<=100)|(t>=149&t<=152)|(t>=201&t<=204)|(t>=253&t<=256)|(t>=305&t<=308)

// Structural dummies
gen DM1 = t>=28 & t<=166
gen DM2 = t>=167 & t<=181
gen DM3 = t>=182 & t<=323
gen DM4 = t>=324 & t<=328

save Data/jecnew,replace

// Summary statistic
sum gr Q L po
