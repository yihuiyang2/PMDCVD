**
** TABLE 1 Basic characteritics
**

sysdir set PLUS “W:\C6_Fang\WMH\Projects\Yihui\IT\Stata”

cd “W:\C6_Fang\WMH\Projects\Yihui\Project 4 PMDCMD\data\”

*append data on atc from population cohort and sibling cohort
use atc, clear
append using atc_sib
save atc_all,replace

* load data on population and sibling cohort
use cohort, clear
drop case id

append using cohort_sib,generate(sib)

rename *, lower

drop fstdt_dm

*merge with data on multi dx
merge m:1 lopnr using multidx_fst

drop if _merge==2
drop _merge

*merge with data on atc
merge 1:1 _setno lopnr sib using atc_all

drop if _merge==2
drop _merge

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
*code missing if before year 2006
replace ht=. if year_match <2006

*use of hormone contraceptives before cohort entry 0-no 1-yes
gen hc=0
replace hc=1 if dt_hc<entry
*code missing if before year 2006
replace hc=. if year_match <2006

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

*recode county
encode lan_match, gen(county)

*code history of multi morbidities
foreach i of num 1/6 {

local w: word `i’ of “dm” “asthma” “menstrual” “pcos” “dyslip” “au”

gen `w’=0
replace `w’=1 if fstdt_`w’<entry
drop fstdt_`w’ `w’_npr

}

*code number of atc
replace count_atc=0 if count_atc==.
tab count_atc

*I classify it into: 0, 1-2, >=3
gen atc=0
replace atc=1 if count_atc==1 | count_atc==2
replace atc=2 if count_atc>=3
*code missing if before year 2007
replace atc=. if year_match <2007


*code BMI: 0-normal 1-underweight 2-overweight 3-obese
gen bmigrp=.
replace bmigrp=0 if bmi >= 18.5 & bmi < 25
replace bmigrp=1 if bmi > 0 & bmi < 18.5
replace bmigrp=2 if bmi >= 25 & bmi < 30
replace bmigrp=3 if bmi >= 30 & bmi < .
*replace missing as 9 among parous women
replace bmigrp=9 if bmi==. & parity > 0


*code smoking: replace missing as 9 among parous women
replace smoke=. if smoke==0
replace smoke=9 if smoke==. & parity > 0

*check interval between delivery year and index year
gen year_childbirth=year(dt_childbirth)
gen interval=year_match-year_childbirth
sum interval if sib==0, detail

sum interval if sib==1, detail

* code exit
gen t1=min(fstdt_cvd1, emigrate_dt, death_dt, td(31dec2022)) if _cc==1
replace t1=min(fstdt_cvd1, fstpmd_dt, emigrate_dt, death_dt, td(31dec2022)) if _cc==0
format t1 %d

*code age at matching
gen age_match=(entry-dt_birth)/365.25

gen age_endfu=(t1-dt_birth)/365.25

* code calendar year at matching
gen ygrp=0 if year_match>=2001 & year_match<=2007
replace ygrp=1 if year_match>=2008 & year_match<=2014
replace ygrp=2 if year_match>=2015 & year_match<=2022


*calculate median and max fu in population and sibling analysis
gen py=(t1-entry)/365.25
replace py=round(py,0.1)

sum py if sib==0, detail
*median: 6.2
*max: 22

sum py if sib==1, detail
*median: 6.3
*max: 22


*calculate median age at entry in population and sibling analysis
sum age_match if sib==0, detail
*median 35.4

sum age_match if sib==1, detail
*median 34.6

*calculate median age at end of follow-up in population and sibling analysis
sum age_endfu if sib==0, detail
*median 42.2

sum age_endfu if sib==1, detail
*median 41.3


* keep variables
keep lopnr _setno _cc _rstime byear county edu marr income parity fland psy ht hc sib age_match age_endfu ygrp year_match bmigrp smoke dm asthma menstrual pcos dyslip au atc

*check missing
misstable sum
*seems fine

* calculate N, PY and IR overall
  matrix A = J(1,15,.)
* pop design, PMDs
count if sib==0 & _cc==1
  *matrix A[1,1] = 0
  matrix A[1,2] = 0
  matrix A[1,4] = r(N)
  scalar N1=r(N)

* pop design, non-PMDs
count if sib==0 & _cc==0
  matrix A[1,7] = r(N)
  scalar N2=r(N)

* sib design, PMDs
count if sib==1 & _cc==1
  matrix A[1,10] = r(N)
  scalar N3=r(N)

* sib design, non-PMDs
count if sib==1 & _cc==0
  matrix A[1,13] = r(N)
  scalar N4=r(N)

* calculate N, PY and IR since 2006 (for HT and HC)
* pop design, PMDs
count if sib==0 & _cc==1 & year_match>=2006
  scalar N5=r(N)

* pop design, non-PMDs
count if sib==0 & _cc==0 & year_match>=2006
  scalar N6=r(N)

* sib design, PMDs
count if sib==1 & _cc==1 & year_match>=2006
  scalar N7=r(N)

* pop design, non-PMDs
count if sib==1 & _cc==0 & year_match>=2006
  scalar N8=r(N)



* calculate N, PY and IR among parous women (for smoking and BMI)
* pop design, PMDs
count if sib==0 & _cc==1 & parity>0
  scalar N9=r(N)

* pop design, non-PMDs
count if sib==0 & _cc==0 & parity>0
  scalar N10=r(N)

* sib design, PMDs
count if sib==1 & _cc==1 & parity>0
  scalar N11=r(N)

* pop design, non-PMDs
count if sib==1 & _cc==0 & parity>0
  scalar N12=r(N)

  * calculate N, PY and IR since 2007 (for accumulated medication)
* pop design, PMDs
count if sib==0 & _cc==1 & year_match>=2007
  scalar N13=r(N)

* pop design, non-PMDs
count if sib==0 & _cc==0 & year_match>=2007
  scalar N14=r(N)

* sib design, PMDs
count if sib==1 & _cc==1 & year_match>=2007
  scalar N15=r(N)

* pop design, non-PMDs
count if sib==1 & _cc==0 & year_match>=2007
  scalar N16=r(N)

* output
matrix X = A

* calculate median and IQR by numeric covars
foreach i of num 1/2 {
local j: word `i’ of “age_match” “age_endfu” “ygrp” “byear” “county” “edu” “marr” “income” “parity” “fland” “psy”  “ht” “hc” “dm” “asthma” “menstrual” “pcos” “dyslip” “au” “atc” “smoke” “bmi”
  matrix A[1,1] = 0
  matrix A[1,2] = `i’

* pop design, PMDs
sum `j’ if sib==0 & _cc==1, detail
  matrix A[1,4] =  r(p50)
  matrix A[1,5] =  r(p25)
  matrix A[1,6] =  r(p75)

* pop design, non-PMDs
sum `j’ if sib==0 & _cc==0, detail
  matrix A[1,7] =  r(p50)
  matrix A[1,8] =  r(p25)
  matrix A[1,9] =  r(p75)
	
* sib design, PMDs
sum `j’ if sib==1 & _cc==1, detail
  matrix A[1,10] =  r(p50)
  matrix A[1,11] =  r(p25)
  matrix A[1,12] =  r(p75)

* sib design, non-PMDs
sum `j’ if sib==1 & _cc==0, detail
  matrix A[1,13] =  r(p50)
  matrix A[1,14] =  r(p25)
  matrix A[1,15] =  r(p75)

* output
matrix X = X\A
}


* calculate N, PY and IR by categorical covars
foreach i of num 3/22 {
local j: word `i’ of  “age_match” “age_endfu” “ygrp” “byear” “county” “edu” “marr” “income” “parity” “fland” “psy”  “ht” “hc” “dm” “asthma” “menstrual” “pcos” “dyslip” “au” “atc” “smoke” “bmi”

* loop levels
levelsof `j’, local(lvl)

foreach k of local lvl {
  matrix A = J(1,15,.)

  matrix A[1,1] = 1
  matrix A[1,2] = `i’
  matrix A[1,3] = `k’

* pop design, PMDs
count if sib==0 & _cc==1  & `j’ == `k’
  matrix A[1,4] = r(N)

* pop design, non-PMDs
count if sib==0 & _cc==0  & `j’ == `k’
  matrix A[1,7] = r(N)
	
* sib design, PMDs
count if sib==1 & _cc==1 & `j’ == `k’
  matrix A[1,10] = r(N)

* sib design, non-PMDs
count if sib==1 & _cc==0 & `j’ == `k’
  matrix A[1,13] = r(N)

* output
matrix X = X\A
}
}

matrix colnames X = cat var lvl c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15

xsvmat X, names(col) saving(../output/table1, replace)

*import data
use ../output/table1,clear

*calculate percentage for categorical variables
replace c5 = c4/N1*100 if cat==1
replace c8 = c7/N2*100 if cat==1
replace c11 = c10/N3*100 if cat==1
replace c14 = c13/N4*100 if cat==1

*calculate percentage for ht and hc
replace c5= c4/N5*100 if cat==1 & (var==12 | var==13)
replace c8= c7/N6*100 if cat==1 & (var==12 | var==13)
replace c11= c10/N7*100 if cat==1 & (var==12 | var==13)
replace c14= c13/N8*100 if cat==1 & (var==12 | var==13)

*calculate percentage for smoke and bmi
replace c5= c4/N9*100 if cat==1 & (var==21 | var==22)
replace c8= c7/N10*100 if cat==1 & (var==21 | var==22)
replace c11= c10/N11*100 if cat==1 & (var==21 | var==22)
replace c14= c13/N12*100 if cat==1 & (var==21 | var==22)

*calculate percentage for medication
replace c5= c4/N13*100 if cat==1 & var==20
replace c8= c7/N14*100 if cat==1 & var==20
replace c11= c10/N15*100 if cat==1 & var==20
replace c14= c13/N16*100 if cat==1 & var==20


*format

foreach i of num 1/4 {

*variable list
local a: word `i’ of “c4” “c7" “c10” “c13"
local b: word `i’ of “c5” “c8" “c11” “c14"
local c: word `i’ of “c6” “c9" “c12” “c15"
local d: word `i’ of “pop_pmd” “pop_nonpmd” “sib_pmd” “sib_nonpmd”

*median(IQR)
gen `d’= “”
replace `d’=string(`a’,“%9.1f”) + ” (” + string(`b’,“%9.1f”) + “-” + string(`c’,“%9.1f”) + “)” if cat==0

*number(%)
replace `d’ = string(`a’,“%9.2gc”) + ” (” + string(`b’,“%9.2f”) + “)” if cat==1

*total number
replace `d’ = string(`a’,“%9.2gc”) if var==0
}

* label
label define varlab 0 “Total” 1 “Age at matching” 2 “Age at end of follow-up” 3 “Calendar year at matching” 4 “Birth year” 5 “County of residence” 6 “Educational level” 7 “Civil status” 8 “Income” 9 “Parity” 10 “Country of birth” 11 “Psychiatric diagnosis” 12 “Use of hormone therapy” 13 “Use of hormone contraceptives” 14 “Diabetes” 15 “Asthma” 16"Menstrual disorders” 17 “PCOS” 18 “Dyslipidemia” 19 “Autoimmune diseases” 20 “Medication” 21 “Smoke” 22 “BMI group”
label values var varlab


*keep variables
keep var lvl pop_pmd pop_nonpmd sib_pmd sib_nonpmd


* output
export excel ../output/table1.xlsx, firstrow(var) replace
