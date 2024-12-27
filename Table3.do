
**
** TABLE 3 PMDs and CVD, by PMD subtypes
**

sysdir set PLUS "W:\C6_Fang\WMH\Projects\Yihui\IT\Stata"
cd "W:\C6_Fang\WMH\Projects\Yihui\Project 4 PMDCMD\data\"


* load data
use cohort, clear
drop case id

append using cohort_sib,generate(sib)

rename *, lower

*code entry
gen entry=_rstime

*birth year
gen byear=year(dt_birth)

gen age=year_match-byear

*cvd dx date
gen fstdt_cvd1=min(fstdt_cvd11,fstdt_cvd12,fstdt_cvd13,fstdt_cvd14,fstdt_cvd21,fstdt_cvd31,fstdt_cvd41,fstdt_cvd42,fstdt_cvd43,fstdt_cvd51,fstdt_cvd52,fstdt_cvd61,fstdt_cvd71)
format fstdt_cvd1 %d


*use of hrt before cohort entry 0-no use 1-yes
gen ht=0
replace ht=1 if dt_ht<entry

*use of hormone contraceptives before cohort entry 0-no 1-yes
gen hc=0
replace hc=1 if dt_hc<entry

*recode birth year to groups
egen byear2=cut(byear),at(1949,1959,1969,1979,1989,1999,2007)
drop byear
rename byear2 byear

*history of psychiatric disorders
gen psy=0
replace psy=1 if dt_psy<_rstime

*recode income to quintiles
egen income2 = cut(income), group(5)
replace income2=9 if income==. 
drop income
rename income2 income

*recode edu
replace edu=9 if edu==.

*recode marr
replace marr=9 if marr==.

encode lan_match, generate(county)

*subtype of PMDs, by cormorbid PND at baseline:0-no PMDs 1-PMDs without PND 2-PMDs with PND
gen pms=0 if _cc==0
replace pms=1 if _cc==1 & dt_pnd >= _rstime
replace pms=2 if _cc==1 & dt_pnd < _rstime

*subtype of PMDs, by age at entry: 0-no PMDs 1-PMDs < 25 2-PMDs >= 25

foreach i of num 25 {

  gen pms`i'=0 if _cc==0
  replace pms`i'=1 if _cc==1 & age < `i'
  replace pms`i'=2 if _cc==1 & age >= `i'

}



* loop 
matrix X = J(1,9,.)

foreach t of num 0/1 {
	foreach k of num 1/2 {
	   foreach j of num 1/2 {

preserve
local sub: word `k' of "pms" "pms25" 

* keep obs
keep if sib==`t' 

keep if `sub'==0 | `sub'==`j'

*when studying subtype by pnd, drop non-parous women at baseline
drop if parity==. & `k'==1

*exclude matching set with only one observation
gen count=1
egen total = total(count), by(_setno)
drop if total==1

*exclude matching set without case
egen total2 = total(_cc), by(_setno)
drop if total2==0


* code exit
gen t1=min(fstdt_cvd1, emigrate_dt, death_dt, td(31dec2022)) if _cc==1
replace t1=min(fstdt_cvd1, fstpmd_dt, emigrate_dt, death_dt, td(31dec2022)) if _cc==0
format t1 %d

* code fail
gen fail=0
replace fail=1 if fstdt_cvd1==t1
replace t1=t1+1 if _rstime==t1

* set time - time from matching
gen id=_n
stset t1, failure(fail==1) enter(time _rstime) exit(time t1) origin(time _rstime) scale(365.25) id(id)


**
** PMD and CVD risk

* code exposure
gen exp=`sub'
replace exp=1 if `sub'==2
	 
matrix A = J(2,9,.)
  
*N and Crude IR
foreach i of num 0/1 { 
  stptime if exp==`i', per(1000)
  matrix A[`i'+1,1]  = r(failures) 
  matrix A[`i'+1,2]  = r(rate)
  matrix A[`i'+1,6] = `i'
  matrix A[`i'+1,7] = `t'
  matrix A[`i'+1,8] = `k'
  matrix A[`i'+1,9] = `j'
}

*Model 3: further adjusted for psychiatric history
stcox i.exp i.byear i.county i.income i.edu i.marr i.fland i.ht i.hc i.psy, strata(_setno) vce(cluster lopnr)
foreach i of num 0/1 { 
  matrix A[`i'+1,3]  = exp(_b[`i'.exp])
  matrix A[`i'+1,4] = exp(_b[`i'.exp] - 1.96 * _se[`i'.exp])
  matrix A[`i'+1,5] = exp(_b[`i'.exp] + 1.96 * _se[`i'.exp])  
}


mat X = X\A	 
restore
}
}
}


*output table
matrix list X
matrix colnames X = n ir b1 lb1 ub1 exp sib sub type
xsvmat X, names(col) saving(../output/table3, replace)


*format
use ../output/table3,clear

*drop the first line
drop if exp==.

*format n/ir
gen N = string(n,"%9.2gc")
gen IR =  string(ir,"%9.2f") 

*format HRs
foreach i of num 1 {

*variable list
local c: word `i' of "HR1"
local r: word `i' of "b1" 
local l: word `i' of "lb1" 
local h: word `i' of "ub1"

*HR
gen `c'     = "1.00"
replace `c' = string(`r',"%9.2f") + " (" + string(`l',"%9.2f") + "-" + string(`h',"%9.2f") + ")" if exp>0
}

*create label for sub
label define sublab 1 "pms" 2 "pms25" 
label values sub sublab

*create label for sib
label define siblab 0 "population analysis" 1 "sibling analysis" 
label values sib siblab


*output
keep sib sub type exp N IR HR*
export excel ../output/table3.xlsx, firstrow(var) replace









