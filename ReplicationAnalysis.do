********************************************************************************
***					Replication of Gebel & Pfeiffer (2010)					 ***
********************************************************************************

clear
clear matrix

version 14
set more off, permanently
set rmsg on, permanently

********************************************************************************
***	Directories																 ***
********************************************************************************

global prog /home/avila/Documents/STUDIUM/03SemMaster/empFiWi/ReplicationHeterogeneusEducReturn/Prog
global do ${prog}/Do
global log ${prog}/Log
global temp ${prog}/Temp
global data ${prog}/Data
global output ${prog}/Output


********************************************************************************
***				 Estimations with different approaches						 ***
********************************************************************************

use ${temp}/sample, clear

********************************************************************************
*** Mincerian type OLS regression (homogenous returns to education)
********************************************************************************

*** Run regression by year
gen ret_OLS=.
gen se_OLS=.
forv n=1984/2006 {
	reg lnw school age age_sq if syear==`n', robust
	replace ret_OLS = _b[school] if syear==`n'
	replace se_OLS = _se[school] if syear==`n'
}
tabstat ret* se_*, by(syear)

*** Regress with interactions of schooling & syear
reg lnw c.school#i.syear age age_sq, robust

* Plot marginal effects by syear
margins syear, dydx(school)
*marginsplot

* Graph estimated return to education per year
sort syear
twoway line ret_OLS syear, ylabel(0 (0.02) 0.14) ytitle("Average return to education") ///
xtitle("") xlabel(1984 (2) 2006)


********************************************************************************
*** CF approach
********************************************************************************

*** Schooling = endogeneous variable

********************************************************************************
*** Regress with OLS per year

gen ret_cf=.
gen se_cf=.

forv n=1984/2006 {

* First stage: Estimate the reduced form of schooling, i.e. regress schooling 
* on all exogeneous variables including the instrument (siblings)
reg school sex age age_sq rural edu_f edu_m occ_f sibl if syear==`n'

* Obtain the residuals
predict v`n', res

* Second stage:  Estimate the structural equation and include the residuals from 
* the reduced form as an additional regressor
reg lnw school sex age age_sq rural edu_f edu_m occ_f v`n' c.v`n'#c.school if syear==`n' // exclude instrument
replace ret_cf = _b[school] if syear==`n'
replace se_cf = _se[school] if syear==`n'

* Endogeneity test: test whether the coefficients on the residuals are statistically significant
test v`n' c.v`n'#c.school														// reject v=0 & v#c.school=0 at 1% significance level for all years, except 1998/1999 at 10% level
}

* Graph estimated return to education per year
sort syear
twoway (line ret_OLS syear) (line ret_cf syear), ylabel(0 (0.02) 0.14) ytitle("Average return to education") ///
xtitle("") xlabel(1984 (2) 2006)

* Marginal effect of schooling
margins, dydx(school) by(syear)

********************************************************************************
*** CMI approach															 ***
********************************************************************************

*** Generalized linear regression model with Poisson distribution

glm school sex age age_sq rural edu_f occ_f edu_m, family(poisson) link(log)

*** Number of observations
local N = e(N)
scalar x = 1/`N'

*** Predict conditional mean on observable characteristics x (expected value of school)
predict cm, mu

*** Predict pearson residuals
predict res_p, pearson

*** Calculate residual as difference between prediction & observed schooling
gen res = school - cm

*** Estimate sigma_sq
egen sigma_sq = mean(res^2)
egen sigma_sq_pearson = mean(res_p^2)

*** Estimate the variance
gen var = sigma_sq*cm
gen var_pears = sigma_sq_pearson*cm

*** Calculate the APE-estimator

* Inner term
gen inner = (res*lnw)/var
gen inner_pears = (res*lnw)/var_pears

* Sum the inner term over all observations
egen sum = total(inner)
egen sum_pears = total(inner_pears)

* Calculate coefficient
gen b_hat = x*sum																// sum is too small (and too many observations)
gen b_hat_pears = x*sum_pears

egen aa = mean((res*lnw)/ (sigma_sq_pearson*cm))

*** Compare
sum school cm res sigma_sq var inner sum
sum b_hat b_hat_pear 

********************************************************************************
*** Repeat per  year

forv n = 1984/2006 {

	*** Generalized linear regression model with Poisson distribution
	di "`n'"
	glm school sex age age_sq rural edu_f occ_f edu_m if syear == `n', family(poisson) link(log)

	*** Number of observations
	local N = e(N)
	scalar x = 1/`N'

	*** Predict conditional mean on observable characteristics x (expected value of school)
	predict cm`n', mu

	*** Predict pearson residuals
	predict res_p`n', pearson

	*** Calculate residual as difference between prediction & observed schooling
	gen res`n' = school - cm`n'

	*** Estimate sigma_sq
	egen sigma_sq`n' = mean(cond(syear==`n' & !mi(res`n'),res`n'^2,.))
	egen sigma_sq_pearson`n' = mean(cond(syear==`n' & !mi(res_p`n'),res_p`n'^2,.))

	*** Estimate the variance
	gen var`n' = sigma_sq`n'*cm`n'
	gen var_pears`n' = sigma_sq_pearson`n'*cm`n'

	*** Calculate the APE-estimator

	* Inner term
	gen inner`n' = (res`n'*lnw)/var`n' if syear==`n'
	gen inner_pears`n' = (res`n'*lnw)/var_pears`n' if syear==`n'

	* Sum the inner term over all observations
	egen sum`n' = total(cond(syear==`n' & !mi(inner`n'),inner`n',.))				// sum is too small (and too many observations)
	egen sum_pears`n' = total(cond(syear==`n'& !mi(inner_pears`n'),inner_pears`n',.))

	* Calculate coefficient
	gen b_hat`n' = x*sum`n'
	gen b_hat_pears`n' = x*sum_pears`n'
}

* Summarize in one  variable
gen ret_cmi = .
forv n = 1984/2006 {
	replace ret_cmi = b_hat_pears`n' if syear==`n'
}

*** Compare APE estimates
sum b_hat1* b_hat2*
sum b_hat_p*

********************************************************************************
*** Bootstrapping
********************************************************************************

********************************************************************************
*** Define Program myCMI

capture program drop myCMI           // drop if alredy created. 
program define myCMI, rclass  
  preserve    // using preserve/restore otherwise variable already defined error 
    bsample   // not sure why or if needed :D 
	
	*** Generalized linear regression model with Poisson distribution
	quietly glm school sex age age_sq rural edu_f occ_f edu_m, family(poisson) link(log)

	*** Predict conditional mean on observable characteristics x 
	* (expected value of school) and pearson residuals
	predict condMean, mu
	predict res_pears, pearson

	*** Calculate residual as difference between prediction & observed schooling
	gen resid = school - condMean

	*** Estimate sigma_sq
	egen sigma_sq_pears = mean(res_pears^2)

	egen bCMI = mean((resid*lnw)/ (sigma_sq_pears*condMean))
	local bCMI=bCMI // doesnt look pretty, but returning bCMI directly returns error. 
	
	return scalar bCMI_return = `bCMI'
  restore
end

********************************************************************************
*** Run Bootstrapping

*** define loop values 
local endYear = 2006    // short period for testing purposes 
local reps = 500          // short reps as well 
gen se_cmi = . 

*** run actual boostrapping 
forv n = 1984/`endYear' {
  bootstrap r(bCMI_return), reps(`reps') seed(42) saving(${temp}/bsCMI_`n', replace): myCMI if syear == `n'
  matrix SE = e(se)
  replace se_cmi = SE[1,1] if syear == `n'
}

*** plot histogram to check results 
capture preserve
  use ${temp}/bsCMI_1990, clear 
  histogram _bs_1
capture restore

********************************************************************************
*** Compare all estimates													 ***
********************************************************************************

*** Graph estimated returns to education per year
sort syear
twoway (line ret_OLS syear) (line ret_cf syear) (line ret_cmi syear), ///
    ylabel(0 (0.02) 0.14) ytitle("Average return to education") ///
    xtitle("") xlabel(1984 (2) 2006) legend(order( 1 "OLS" 2 "APE(CF)" 3 "APE(CMI)"))

graph export $output/graph_OLS_CMI_CF.png, replace

**** Table
tabstat *OLS* *cf* *cmi*, by(syear) format(%9.4f) save


