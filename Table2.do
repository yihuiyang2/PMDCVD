
**
** TABLE 2 PMDs and CVD
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


* loop 
matrix X = J(1,13,.)

foreach t of num 0/1 {

preserve
* keep obs
keep if sib==`t' 

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
gen exp=_cc
	 
matrix A = J(2,13,.)
  
*N and Crude IR
foreach i of num 0/1 { 
  stptime if exp==`i', per(1000)
  matrix A[`i'+1,1]  = r(failures) 
  matrix A[`i'+1,2]  = r(rate)
  matrix A[`i'+1,12] = `i'
  matrix A[`i'+1,13] = `t'
}

*Model 1: adjusted for birth year and county
stcox i.exp i.byear i.county, strata(_setno) vce(cluster lopnr)
foreach i of num 0/1 { 
  matrix A[`i'+1,3] = exp(_b[`i'.exp])
  matrix A[`i'+1,4] = exp(_b[`i'.exp] - 1.96 * _se[`i'.exp])
  matrix A[`i'+1,5] = exp(_b[`i'.exp] + 1.96 * _se[`i'.exp])  
}
  
*Model 2: further adjusted for use of HRT, use of HC, educational level, country of birth, civil status, and income
stcox i.exp i.byear i.county i.income i.edu i.marr i.fland i.ht i.hc, strata(_setno) vce(cluster lopnr)
foreach i of num 0/1 { 
  matrix A[`i'+1,6] = exp(_b[`i'.exp])
  matrix A[`i'+1,7] = exp(_b[`i'.exp] - 1.96 * _se[`i'.exp])
  matrix A[`i'+1,8] = exp(_b[`i'.exp] + 1.96 * _se[`i'.exp])  
}
  
*Model 3: further adjusted for psychiatric history
stcox i.exp i.byear i.county i.income i.edu i.marr i.fland i.ht i.hc i.psy, strata(_setno) vce(cluster lopnr)
foreach i of num 0/1 { 
  matrix A[`i'+1,9]  = exp(_b[`i'.exp])
  matrix A[`i'+1,10] = exp(_b[`i'.exp] - 1.96 * _se[`i'.exp])
  matrix A[`i'+1,11] = exp(_b[`i'.exp] + 1.96 * _se[`i'.exp])  
}


mat X = X\A	 
restore
}


*output table
matrix list X
matrix colnames X = n ir b1 lb1 ub1 b2 lb2 ub2 b3 lb3 ub3 exp sib
xsvmat X, names(col) saving(../output/table2, replace)


*format
use ../output/table2,clear

*drop the first line
drop if exp==.

*format n/ir
gen N = string(n,"%9.2gc")
gen IR =  string(ir,"%9.2f") 

*format HRs
foreach i of num 1/3 {

*variable list
local c: word `i' of "HR1" "HR2" "HR3" 
local r: word `i' of "b1" "b2" "b3"
local l: word `i' of "lb1" "lb2" "lb3"
local h: word `i' of "ub1" "ub2" "ub3" 

*HR
gen `c'     = "1.00"
replace `c' = string(`r',"%9.2f") + " (" + string(`l',"%9.2f") + "-" + string(`h',"%9.2f") + ")" if exp>0
}

*create label for sib
label define siblab 0 "population analysis" 1 "sibling analysis" 
label values sib siblab

*output
keep sib exp N IR HR*
export excel ../output/table2.xlsx, firstrow(var) replace









