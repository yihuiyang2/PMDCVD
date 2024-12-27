
**
** TABLE 4 PMDs and CVD, by subtypes of CVD
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

*first hypertensive disease dx date
gen fstdt_cvd50=min(fstdt_cvd51,fstdt_cvd52)
format fstdt_cvd50 %d

*first Cerebrovascular disease dx date
gen fstdt_cvd10=min(fstdt_cvd11,fstdt_cvd12,fstdt_cvd13,fstdt_cvd14)
format fstdt_cvd10 %d

*put cardiac arrest to other types of CVD
replace fstdt_cvd71=min(fstdt_cvd71,fstdt_cvd43)

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

* recode parity
recode parity (1=1) (2=2) (3=3) (4/17=4)
replace parity=0 if parity==.

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
matrix X = J(1,8,.)
foreach j of num 1/14 {
		foreach t of num 0/1 {
local s: word `j' of "fstdt_cvd50" "fstdt_cvd51" "fstdt_cvd52" "fstdt_cvd61" "fstdt_cvd10" "fstdt_cvd11" "fstdt_cvd12" "fstdt_cvd13" "fstdt_cvd14" "fstdt_cvd21" "fstdt_cvd31"  "fstdt_cvd41" "fstdt_cvd42" "fstdt_cvd71"    

preserve
* keep obs
keep if sib==`t' 

* code exit
gen t1=min(`s', emigrate_dt, death_dt, td(31dec2022)) if _cc==1
replace t1=min(`s', fstpmd_dt, emigrate_dt, death_dt, td(31dec2022)) if _cc==0
format t1 %d

* code fail
gen fail=0
replace fail=1 if `s'==t1
replace t1=t1+1 if _rstime==t1

* set time - time from matching
gen id=_n
stset t1, failure(fail==1) enter(time _rstime) exit(time t1) origin(time _rstime) scale(365.25) id(id)


**
** PMD and CVD risk

* code exposure
gen exp=_cc
	 
matrix A = J(2,8,.)
  
*N and Crude IR
foreach i of num 0/1 { 
  stptime if exp==`i', per(1000)
  matrix A[`i'+1,1]  = r(failures) 
  matrix A[`i'+1,2]  = r(rate)
  matrix A[`i'+1,6] = `i'
  matrix A[`i'+1,7] = `j'
  matrix A[`i'+1,8] = `t'
}

  
*Model 3
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

*output table
matrix list X
matrix colnames X = n ir b3 lb3 ub3 exp outcome sib
xsvmat X, names(col) saving(../output/table4, replace)

*format
use ../output/table4,clear

*drop the first line
drop if exp==.

*format n/ir
gen N = string(n,"%9.2gc")
gen IR =  string(ir,"%9.2f") 

*HR
gen HR3 = "1.00"
replace HR3 = string(b3,"%9.2f") + " (" + string(lb3,"%9.2f") + "-" + string(ub3,"%9.2f") + ")" if exp>0


*create label for outcome
label define outlab 1 "Hypertensive diseases" 2 "Essential hypertension" 3 "Other hypertensive disease" 4 "Ischemic heart disease"  5 "Cerebrovascular disease" 6 "Subarachnoid hemorrhage" 7 "Hemorrhagic stroke" 8 "Ischemic stroke" 9 "Other cerebrovascular disease" 10  "Emboli and thrombosis" 11 "Heart failure" 12 "Arrhythmia" 13 "Conduction disorder" 14 "Other type of CVD"
label values outcome outlab

*create label for sib
label define siblab 0 "population analysis" 1 "sibling analysis" 
label values sib siblab


*output
keep sib exp outcome N IR HR*
export excel ../output/table4.xlsx, firstrow(var) replace








