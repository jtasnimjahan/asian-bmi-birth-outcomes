/******************************************************************************
*** 01_data_prep.sas
***
*** Project : Hidden overweight and adverse birth outcomes among Asian mothers
***           in the United States, 2016-2023
*** Purpose : (A) Read the NVSS Natality public-use files, 2016-2023, pool
***               them (ana.allyears) and build the analytic file
***               ana.analytic_pooled.
***           (B) Derive the analysis variables used by 02 and 03 and save
***               the Asian + White analysis file ana.wa.
***           (C) Compare excluded and included records.
***           (D) Table 1 with absolute standardized differences.
*** Produces: Figure 1 counts (Part A); Table 1 (Part D);
***           eTable 2 (Part C)
*** Software: SAS 9.4
*** Author  : Jabin Tasnim Jahan, University at Albany, SUNY
*** Contact : jtasnimjahan@albany.edu
*** Version : 1.0 (November 2026), code for the submitted manuscript
******************************************************************************/

options nocenter ls=200 ps=max;

/* Set these paths before running */
%let RAWDIR = /path/to/nvss/files;     /* folder with the raw NVSS .txt files */
libname ana "/path/to/your/data";      /* folder for SAS datasets             */
%let OUTDIR = /path/to/your/output;    /* folder for results                  */


/*============================================================================
  PART A. READ, POOL AND BUILD ana.analytic_pooled
============================================================================*/

/* NVSS public-use file names, as downloaded from CDC */
%let f2016 = &RAWDIR./Nat2016PublicUS.c20170517.r20190620.txt;
%let f2017 = &RAWDIR./Nat2017PublicUS.c20180516.r20180808.txt;
%let f2018 = &RAWDIR./Nat2018PublicUS.c20190509.r20190717.txt;
%let f2019 = &RAWDIR./Nat2019PublicUS.c20200506.r20200915.txt;
%let f2020 = &RAWDIR./Nat2020PublicUS.c20210506.r20210812.txt;
%let f2021 = &RAWDIR./Nat2021US.txt;
%let f2022 = &RAWDIR./Nat2022PublicUS.c20230504.r20230822.txt;
%let f2023 = &RAWDIR./Nat2023PublicUS.c20240509.r20240724.txt;

/*---------------------------------------------------------------------------
  A1. Read each year. Field positions are the same for 2016-2023.
      BMI is read with a 4. informat because the raw field contains an
      explicit decimal point in every year.
---------------------------------------------------------------------------*/
%macro readnat(yr, filepath);
  data yr&yr;
    infile "&filepath" lrecl=2000 missover;
    input
      @75   mager        2.    /* maternal age, single year           */
      @84   mbstate_rec  1.    /* nativity: 1 US, 2 foreign, 3 unknown */
      @108  mrace15      2.    /* race recode 15 (04 = Asian Indian)   */
      @117  mracehisp    1.    /* race/Hispanic origin                 */
      @120  dmar         1.    /* marital status                       */
      @124  meduc        1.    /* maternal education                   */
      @179  lbo_rec      1.    /* live birth order recode              */
      @227  precare5     1.    /* trimester prenatal care began        */
      @269  cig_rec     $1.    /* smoking recode Y/N/U                 */
      @283  bmi          4.    /* prepregnancy BMI                     */
      @287  bmi_r        1.    /* BMI recode (9 = unknown)             */
      @316  rf_ghype    $1.    /* gestational hypertension             */
      @317  rf_ehype    $1.    /* eclampsia                            */
      @407  rdmeth_rec   1.    /* delivery method recode               */
      @436  pay_rec      1.    /* payment source recode                */
      @454  dplural      1.    /* plurality (1 = singleton)            */
      @475  sex         $1.    /* infant sex                           */
      @499  oegest_comb  2.    /* obstetric estimate of gestation      */
      @504  dbwt         4.    /* birthweight, grams                   */
    ;
    year = &yr;
  run;

  proc append base=ana.allyears data=yr&yr force; run;
  proc datasets lib=work nolist; delete yr&yr; quit;
%mend readnat;

proc datasets lib=ana nolist; delete allyears; quit;

%readnat(2016, &f2016);
%readnat(2017, &f2017);
%readnat(2018, &f2018);
%readnat(2019, &f2019);
%readnat(2020, &f2020);
%readnat(2021, &f2021);
%readnat(2022, &f2022);
%readnat(2023, &f2023);

title "Check: records by year";
proc freq data=ana.allyears;
  tables year / nocum;
run;

/*---------------------------------------------------------------------------
  A2. Exclusions, applied in order (Figure 1)
---------------------------------------------------------------------------*/
data step1;                                   /* singleton births */
  set ana.allyears;
  if dplural = 1;
run;

data step2;                                   /* valid prepregnancy BMI */
  set step1;
  if 13.0 <= bmi <= 69.9 and bmi_r ne 9;
run;

data analytic_pooled;                         /* valid age, race/ethnicity, nativity */
  set step2;
  if 12 <= mager <= 50;
  if mracehisp in (1,2,3,4,5,6,7);
  if mbstate_rec in (1,2);
run;

proc sql;
  title "Figure 1: exclusion flow";
  select "All records 2016-2023"           as step length=40, count(*) as n format=comma12. from ana.allyears
  union all
  select "Singleton births",                  count(*) from step1
  union all
  select "Valid prepregnancy BMI",            count(*) from step2
  union all
  select "Eligible singleton births",         count(*) from analytic_pooled;
quit;
title;

/*---------------------------------------------------------------------------
  A3. Analysis variables
---------------------------------------------------------------------------*/
data ana.analytic_pooled;
  set analytic_pooled;

  /* Maternal age (ref 25-29) */
  length agecat $5;
  if      mager < 20 then agecat = '<20';
  else if mager <= 24 then agecat = '20-24';
  else if mager <= 29 then agecat = '25-29';
  else if mager <= 34 then agecat = '30-34';
  else if mager <= 39 then agecat = '35-39';
  else                     agecat = '>=40';

  /* Standard WHO BMI categories (ref normal), used in Table 2 */
  length bmi_std $7;
  if      bmi < 18.5 then bmi_std = 'under';
  else if bmi < 25.0 then bmi_std = 'normal';
  else if bmi < 30.0 then bmi_std = 'over';
  else if bmi < 35.0 then bmi_std = 'obese1';
  else if bmi < 40.0 then bmi_std = 'obese2';
  else                    bmi_std = 'obese3';

  /* Race/ethnicity */
  length raceeth $8;
  select (mracehisp);
    when (1) raceeth = 'NHWhite';
    when (2) raceeth = 'NHBlack';
    when (3) raceeth = 'NHAIAN';
    when (4) raceeth = 'NHAsian';
    when (5) raceeth = 'NHNHOPI';
    when (7) raceeth = 'Hispanic';
    otherwise raceeth = 'Other';
  end;
  asian = (raceeth = 'NHAsian');

  /* Nativity */
  length nativity $11;
  if      mbstate_rec = 1 then nativity = 'USborn';
  else if mbstate_rec = 2 then nativity = 'Foreignborn';

  /* Outcomes */
  if oegest_comb in (.,99) then ptb = .;
  else ptb = (oegest_comb < 37);

  if dbwt in (.,9999) then lbw = .;
  else lbw = (dbwt < 2500);

  if      rdmeth_rec in (3,4) then cesarean = 1;
  else if rdmeth_rec in (1,2) then cesarean = 0;
  else cesarean = .;

  if upcase(rf_ghype) = 'Y' or upcase(rf_ehype) = 'Y' then htn_disorder = 1;
  else if upcase(rf_ghype) = 'N' and upcase(rf_ehype) = 'N' then htn_disorder = 0;
  else htn_disorder = .;

  /* Smoking */
  if      upcase(cig_rec) = 'Y' then smoke = 1;
  else if upcase(cig_rec) = 'N' then smoke = 0;
  else smoke = .;
run;

proc datasets lib=work nolist; delete step1 step2 analytic_pooled; quit;

title "Check: mean BMI by year (all years should be about 26-29)";
proc means data=ana.analytic_pooled n mean min max maxdec=2;
  class year;
  var bmi;
run;

title "Check: analytic sample by race/ethnicity";
proc freq data=ana.analytic_pooled;
  tables raceeth / nocum;
run;
title;


/*============================================================================
  PART B. ANALYSIS FILE: NH Asian + NH White mothers with analysis variables
============================================================================*/

data ana.wa;
  set ana.analytic_pooled;
  where raceeth in ('NHWhite','NHAsian') and 13.0 <= bmi <= 69.9;

  /* Asian-specific BMI categories (WHO action points 23.0 and 27.5) */
  length bmi5 $7;
  if      bmi < 18.5 then bmi5 = 'under';
  else if bmi < 23.0 then bmi5 = 'normal';
  else if bmi < 25.0 then bmi5 = 'hidden';
  else if bmi < 27.5 then bmi5 = 'visover';
  else                    bmi5 = 'obese';
  bmi5_under   = (bmi5 = 'under');
  bmi5_hidden  = (bmi5 = 'hidden');
  bmi5_visover = (bmi5 = 'visover');
  bmi5_obese   = (bmi5 = 'obese');

  /* Slope-change term at BMI 23.0 (eTable 4) */
  hinge23 = max(bmi - 23, 0);

  /* Missing-indicator marital status (eTable 6) */
  dmar_mi = dmar;
  if missing(dmar) then dmar_mi = 9;

  /* Hypertensive disorders excluding eclampsia (eTable 6) */
  if      upcase(rf_ghype) = 'Y' then htn_noecl = 1;
  else if upcase(rf_ghype) = 'N' then htn_noecl = 0;
  else htn_noecl = .;

  /* Foreign-born indicator (three-way interaction, eTable 5) */
  fb = (nativity = 'Foreignborn');

  /* Detailed Asian subgroup */
  length asiangrp $12;
  if asian = 1 then select (mrace15);
    when (4)  asiangrp = 'AsianIndian';
    when (5)  asiangrp = 'Chinese';
    when (6)  asiangrp = 'Filipino';
    when (7)  asiangrp = 'Japanese';
    when (8)  asiangrp = 'Korean';
    when (9)  asiangrp = 'Vietnamese';
    when (10) asiangrp = 'OtherAsian';
    otherwise asiangrp = 'Check';
  end;
run;

title "Check: Asian subgroup coding (no 'Check' rows expected)";
proc freq data=ana.wa;
  where asian = 1;
  tables asiangrp / missing nocum;
run;
title;


/*============================================================================
  PART C. eTable 2: EXCLUDED VERSUS INCLUDED RECORDS
  Uses all 2016-2023 records (ana.allyears) with the three exclusion steps
  combined into one flag.
============================================================================*/

data allflag;
  set ana.allyears;
  included = 1;
  if dplural ne 1 then included = 0;                              /* multiple gestation */
  if not (13.0 <= bmi <= 69.9) or bmi_r = 9 then included = 0;    /* missing/implausible BMI */
  if not (12 <= mager <= 50) then included = 0;                   /* missing maternal age */
  if mracehisp not in (1,2,3,4,5,6,7) then included = 0;          /* unknown race/ethnicity */
  if mbstate_rec not in (1,2) then included = 0;                  /* unknown nativity */

  if dbwt in (.,9999) then bwt_lt2500 = .;
  else bwt_lt2500 = (dbwt < 2500);

  if      rdmeth_rec in (3,4) then ces_raw = 1;
  else if rdmeth_rec in (1,2) then ces_raw = 0;
  else ces_raw = .;
run;

title "eTable 2: records by inclusion";
proc freq data=allflag;
  tables included / missing nocum;
run;

title "eTable 2: maternal age and birthweight by inclusion";
proc means data=allflag mean std maxdec=1;
  class included;
  var mager dbwt;
run;

title "eTable 2: caesarean delivery and low birthweight by inclusion";
proc freq data=allflag;
  tables included*(ces_raw bwt_lt2500) / nocol nopercent missing;
run;

title "eTable 2: race/Hispanic origin by inclusion";
proc freq data=allflag;
  tables included*mracehisp / nocol nopercent missing;
run;
title;


/*============================================================================
  PART D. TABLE 1: CHARACTERISTICS OF NH ASIAN AND NH WHITE MOTHERS
  n (%) by race, and |d| = absolute standardized difference:
    categories:  |p1 - p2| / sqrt((p1(1-p1) + p2(1-p2)) / 2)
    means:       |m1 - m2| / sqrt((s1**2 + s2**2) / 2)
  Unknown and missing values are shown as their own rows. All category
  variables are character, so the macro can join Asian and White rows.
============================================================================*/

data t1;
  set ana.wa;
  length educ5 $22 lbo4 $10 precare $16 pay $8 marital $9 smokec $7;
  if      meduc in (1,2) then educ5 = '1 Less than HS';
  else if meduc = 3      then educ5 = '2 HS/GED';
  else if meduc in (4,5) then educ5 = '3 Some college';
  else if meduc = 6      then educ5 = '4 Bachelor';
  else if meduc in (7,8) then educ5 = '5 Graduate';
  else                        educ5 = '9 Unknown';
  if      lbo_rec in (1,2,3)     then lbo4 = put(lbo_rec, 1.);
  else if 4 <= lbo_rec <= 8      then lbo4 = '4+';
  else                                lbo4 = '9 Unknown';
  select (precare5);
    when (1) precare = '1 1st trimester';
    when (2) precare = '2 2nd trimester';
    when (3) precare = '3 3rd trimester';
    when (4) precare = '4 None';
    otherwise precare = '9 Unknown';
  end;
  select (pay_rec);
    when (1) pay = 'Medicaid';
    when (2) pay = 'Private';
    when (3) pay = 'Self-pay';
    when (4) pay = 'Other';
    otherwise pay = 'Unknown';
  end;
  if      dmar = 1 then marital = 'Married';
  else if dmar = 2 then marital = 'Unmarried';
  else                  marital = 'Missing';
  if      smoke = 1 then smokec = 'Yes';
  else if smoke = 0 then smokec = 'No';
  else                   smokec = 'Unknown';
run;

/* Categorical rows */
proc datasets lib=work nolist; delete table1; quit;

%macro smdcat(var);
  proc freq data=t1 noprint;
    tables raceeth*&var / missing outpct out=_f(keep=raceeth &var count pct_row);
  run;
  proc sql;
    create table _s as
    select "&var" as variable length=12,
           coalescec(a.&var, w.&var) as level length=22,
           a.count as n_asian, a.pct_row as pct_asian,
           w.count as n_white, w.pct_row as pct_white
    from _f(where=(raceeth='NHAsian')) as a
         full join _f(where=(raceeth='NHWhite')) as w
         on a.&var = w.&var;
  quit;
  data _s;
    set _s;
    pa = pct_asian / 100; pw = pct_white / 100;
    abs_d = abs(pa - pw) / sqrt((pa*(1-pa) + pw*(1-pw)) / 2);
    drop pa pw;
  run;
  proc append base=table1 data=_s force; run;
%mend smdcat;

%smdcat(agecat);
%smdcat(nativity);
%smdcat(educ5);
%smdcat(marital);
%smdcat(lbo4);
%smdcat(precare);
%smdcat(pay);
%smdcat(smokec);

title "Table 1: categorical characteristics, n (%) and |d|";
proc print data=table1 noobs;
  format n_asian n_white comma12. pct_asian pct_white 5.1 abs_d 5.2;
run;

/* Continuous rows: maternal age and BMI */
proc means data=t1 noprint nway;
  class raceeth;
  var mager bmi;
  output out=_m mean=m_age m_bmi std=s_age s_bmi;
run;

data table1_means;
  merge _m(where=(raceeth='NHAsian') rename=(m_age=ma_age m_bmi=ma_bmi s_age=sa_age s_bmi=sa_bmi))
        _m(where=(raceeth='NHWhite') rename=(m_age=mw_age m_bmi=mw_bmi s_age=sw_age s_bmi=sw_bmi));
  d_age = abs(ma_age - mw_age) / sqrt((sa_age**2 + sw_age**2) / 2);
  d_bmi = abs(ma_bmi - mw_bmi) / sqrt((sa_bmi**2 + sw_bmi**2) / 2);
  keep ma_age sa_age mw_age sw_age d_age ma_bmi sa_bmi mw_bmi sw_bmi d_bmi;
run;

title "Table 1: mean (SD) maternal age and BMI, Asian vs White, and |d|";
proc print data=table1_means noobs;
  format ma_: sa_: mw_: sw_: 5.1 d_: 5.2;
run;
title;
