*recreate 2012 results for Entero qPCR code to be used in R markdown
*read in datafile 
set more off

import delimited "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/neearcoreall.csv", varnames(1) clear


drop if venfest==1
keep if gibase~=1 & vomitbase~=1


gen exp=log10mean_epcr
replace exp=0 if bodycontact==0
gen swimtemp=bodycontact

preserve

xi : glm hcgi swimtemp exp i.beach milescat age, link(identity) family(binomial) vce(cluster hh)
est store all

predictnl ar =_b[swimtemp]*swimtemp+_b[exp]*exp if e(sample), ci(arlo arup) se(arse)
keep if swimtemp==1
qui bysort exp (ar): keep if _n==1
drop if missing(ar)

rename exp count
keep ar arlo arup arse count

gen gmean=10^count


foreach var of varlist ar arlo arup {
	replace `var'=`var'*1000
	}
	
	
*figure- all subjects

set scheme s2manual
#delimit ;
line ar arlo arup gmean, pstyle(p1 p2 p2) sort clwidth(medthick medthick medthick) xscale(log)
clpattern(solid dash dash) graphregion(fcolor(white))
xlabel(10 100 1000, grid) xmtick(20 30 40 50 60 70 80 90 200 300 400 500 600 700 800 900 2000 3000 4000 5000, grid)
xtitle("{it: Enterococcus} CCE/100ml (ddCT)" "daily geometric mean", size(small))
ytitle("Swimming-associated GI illness (x 1000)", size(small)) ylabel(#15, angle(45) grid) ymtick(##3)
legend(ring(0) pos(11) symxsize(4) forcesize  cols(1) size(small) order(1 2)
lab(1 "Swimming-associated illness") lab(2 "95% Confidence bound"));
#delimit cr

restore

preserve
xi : glm hcgi swimtemp exp i.beach milescat meanbathers if age<=10, link(identity) family(binomial) vce(cluster hh)
est store under10


predictnl ar =_b[swimtemp]*swimtemp+_b[exp]*exp if e(sample), ci(arlo arup) se(arse)
keep if swimtemp==1
qui bysort exp (ar): keep if _n==1
drop if missing(ar)



rename exp count
keep ar arlo arup arse count

gen gmean=10^count


foreach var of varlist ar arlo arup {
	replace `var'=`var'*1000
	}
	
*figure- 10 and under

set scheme s2manual
#delimit ;
line ar arlo arup gmean, pstyle(p1 p2 p2) sort clwidth(medthick medthick medthick) xscale(log)
clpattern(solid dash dash) graphregion(fcolor(white))
xlabel(10 100 1000, grid) xmtick(20 30 40 50 60 70 80 90 200 300 400 500 600 700 800 900 2000 3000 4000 5000, grid)
xtitle("{it: Enterococcus} CCE/100ml (ddCT)" "daily geometric mean", size(small))
ytitle("Swimming-associated GI illness (x 1000)", size(small)) ylabel(#15, angle(45) grid) ymtick(##3)
legend(ring(0) pos(11) symxsize(4) forcesize  cols(1) size(small) order(1 2)
lab(1 "Swimming-associated illness") lab(2 "95% Confidence bound"));
# delimit cr

