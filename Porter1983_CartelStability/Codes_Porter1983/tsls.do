clear
cd /Users/apple/Dropbox/Documents/Resources/06_Industrial_Organization/04_CU_XiaoEmprical_IO/Assignment1/solution/Codes

use Data/jecnew

// Demand side
ivregress 2sls lnQ L (lngr=m1-m12)
// Supply side
ivregress 2sls lngr L DM1-DM4 po (lnQ=m1-m12)
