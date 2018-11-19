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
*** Back to Analysis
