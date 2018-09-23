********************************************************************************
* Overview of US Economy
********************************************************************************
clear
set more off
set scheme s2mono, permanently
local files : dir "C:\Users\13695\Documents\MATLAB\Macro\hw3_data" files "*.csv"

cd "C:\Users\13695\Documents\MATLAB\Macro\hw3_data"

foreach file in `files' {
	import delimited `file', clear
	split date, p(-) gen(d)
	destring d1 d2 d3, replace
	gen edate = mdy(d2, d3, d1)
	format edate %td
	local name= substr("`file'",1,strpos("`file'", ".")-1)
	save "`name'.dta", replace
	/* do your stuff here */
}

* Q1-3: Different Fractions of GDP

use GDP.dta, clear
tsset edate
merge 1:1 edate using PCE.dta, keep(3) nogen
merge 1:1 edate using GPDI.dta, keep(3) nogen
merge 1:1 edate using GTE.dta, keep(3) nogen
rename w068rcq027sbea gte
merge 1:1 edate using PINCOME.dta, keep(3) nogen
gen cons_gdp=pce/gdp
gen inv_gdp=gpdi/gdp
gen exp_gdp=gte/gdp
gen labor_gdp=pincome/gdp
gen cap_gdp=1-labor_gdp
summarize *_*

tsline cons_gdp inv_gdp exp_gdp
graph export "Figure/gdp_spending_approach.png", width(1200) replace
tsline labor_gdp cap_gdp
graph export "Figure/gdp_income_approach.png", width(1200) replace

keep *_*
unab a: _all
foreach k in `a'{
	gen `k'1=`k'
	gen `k'2=`k'
	gen `k'3=`k'
	gen `k'4=`k'
}
collapse (mean)(*1) (sem)(*2) (max)(*3) (min)(*4)
gen id=_n
reshape long `a', i(id) j(n)
tostring id, replace
replace id="mean" if n==1
replace id="se" if n==2
replace id="max" if n==3
replace id="min" if n==4
export delimited id *_* using "Table/gdp_share", replace

* Q4-5: Labor share and capital share

use LABSHARE.dta, clear
tsset edate
rename labshpusa156nrug labor
gen cap=1-labor
tsline labor cap
graph export "Figure/gdp_income_approach.png", width(1200) replace

keep labor cap
unab a: _all
foreach k in `a'{
	gen `k'1=`k'
	gen `k'2=`k'
	gen `k'3=`k'
	gen `k'4=`k'
}
collapse (mean)(*1) (sem)(*2) (max)(*3) (min)(*4)
gen id=_n
reshape long `a', i(id) j(n)
tostring id, replace
replace id="mean" if n==1
replace id="se" if n==2
replace id="max" if n==3
replace id="min" if n==4
export delimited id lab cap using "Table/gdp_lk", replace


* Q6-7: Growth rate of gdp and consumption

use GDP.dta, clear
merge 1:1 edate using POP.dta, keep(3) nogen
merge 1:1 edate using PCE.dta, keep(3) nogen
gen n=_n
tsset n
rename b230rc0q173sbea pop
gen y_capita=gdp/pop
gen c_capita=pce/pop
gen d_gdp=d.y_capita
gen r_gdp=d_gdp/ y_capita
gen d_cons=d.c_capita
gen r_cons=d_cons/ c_capita
tsset edate
tsline r_gdp r_cons
graph export  "Figure/growth_of_gdp_and_consumption.png", width(1200) replace
keep *_*
local a "r_gdp r_cons"
foreach k in `a'{
	gen `k'1=`k'
	gen `k'2=`k'
	gen `k'3=`k'
	gen `k'4=`k'
}
collapse (mean)(*1) (sem)(*2) (max)(*3) (min)(*4)
gen id=_n
reshape long `a', i(id) j(n)
tostring id, replace
replace id="mean" if n==1
replace id="se" if n==2
replace id="max" if n==3
replace id="min" if n==4
export delimited id *_* using "Table/gdp_contribution", replace


* Q8-9: Unemployment

use UNRATE.dta, clear
tsset edate
tsline unrate
graph export  "Figure/unemployment_rate.png", width(1200) replace

use UEMPMEAN.dta, clear
merge 1:1 edate using UEMPMED.dta, keep(3) nogen
tsset edate
tsline uempmean uempmed
graph export "Figure/unemployment_rate.png", width(1200) replace

local a "uempmean uempmed"
foreach k in `a'{
	gen `k'1=`k'
	gen `k'2=`k'
	gen `k'3=`k'
	gen `k'4=`k'
}
collapse (mean)(*1) (sem)(*2) (max)(*3) (min)(*4)
gen id=_n
reshape long `a', i(id) j(n)
tostring id, replace
replace id="mean" if id=="1"
replace id="se" if id=="2"
replace id="max" if id=="3"
replace id="min" if id=="4"
export delimited id `a' using "Table/unemployment", replace



