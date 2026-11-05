/******************************************************************************
*** 03_sensitivity.sas
***
*** Project : Hidden overweight and adverse birth outcomes among Asian mothers
***           in the United States, 2016-2023
*** Purpose : Sensitivity and robustness analyses
*** Input   : ana.wa (from 01_data_prep.sas), ana.analytic_pooled
*** Produces: Section 1  eTable 6 Panel A (missing-indicator and
***                      marital-omitted models); eTable 3 Panel B Ns
***           Section 2  eTable 6 Panel A (multiple imputation, HTN)
***           Section 3  eTable 6 Panel B (eclampsia excluded; 2022-2023)
***           Section 4  eTable 5 Panel A (race x BMI spline tests by
***                      nativity, excluding 2020-2021, three-way test)
***           Section 5  eTable 6 Panel C (E-values)
*** Software: SAS 9.4
*** Author  : Jabin Tasnim Jahan, University at Albany, SUNY
*** Contact : jtasnimjahan@albany.edu
*** Version : 1.0 (November 2026), code for the submitted manuscript
******************************************************************************/

options nocenter ls=200 ps=max mprint;

/* Set these two paths before running */
libname ana "/path/to/your/data";
%let OUTDIR = /path/to/your/output;

ods rtf file="&OUTDIR./03_sensitivity_results.rtf" style=journal;

%let COV   = agecat meduc dmar lbo_rec precare5 pay_rec sex smoke nativity year;
%let COVNM = agecat meduc lbo_rec precare5 pay_rec sex smoke nativity year;   /* no marital */
%let CLSB  = agecat(ref='25-29') meduc(ref='3') lbo_rec(ref='1') precare5(ref='1')
             pay_rec(ref='2') sex(ref='F') smoke(ref='0') nativity(ref='USborn');
%let CLS   = &CLSB dmar(ref='1') year(ref='2016');
%let BMIX  = bmi5_under bmi5_hidden bmi5_visover bmi5_obese;
/* Nativity-stratified models: nativity is constant, so it is dropped */
%let COVNN = agecat meduc dmar lbo_rec precare5 pay_rec sex smoke year;
%let CLSNN = agecat(ref='25-29') meduc(ref='3') lbo_rec(ref='1') precare5(ref='1')
             pay_rec(ref='2') sex(ref='F') smoke(ref='0') dmar(ref='1') year(ref='2016');

proc datasets lib=work nolist kill; quit;

data wa;     set ana.wa;              run;
data asian5; set wa; where asian = 1; run;

%macro bmi5est;
  estimate "Underweight <18.5 vs 18.5-22.9"    bmi5_under   1 / exp cl;
  estimate "Hidden 23.0-24.9 vs 18.5-22.9"     bmi5_hidden  1 / exp cl;
  estimate "Visible OW 25.0-27.4 vs 18.5-22.9" bmi5_visover 1 / exp cl;
  estimate "Obese >=27.5 vs 18.5-22.9"         bmi5_obese   1 / exp cl;
%mend bmi5est;


/*============================================================================
  SECTION 1. MISSING MARITAL STATUS: MISSING INDICATOR AND OMITTED
============================================================================*/

%macro missdata(outcome=, label=);
  title "Section 1a: missing-indicator marital status - &label";
  proc glimmix data=asian5 empirical=classical;
    class &CLSB dmar_mi(ref='1') year(ref='2016') asiangrp(ref='AsianIndian');
    model &outcome = &BMIX asiangrp agecat meduc dmar_mi lbo_rec precare5
                     pay_rec sex smoke nativity year / dist=poisson link=log;
    random _residual_;
    %bmi5est;
  run;

  title "Section 1b: marital status omitted - &label";
  proc glimmix data=asian5 empirical=classical;
    class &CLSB year(ref='2016') asiangrp(ref='AsianIndian');
    model &outcome = &BMIX asiangrp &COVNM / dist=poisson link=log;
    random _residual_;
    %bmi5est;
  run;
%mend missdata;

%missdata(outcome=htn_disorder, label=Hypertensive disorders);
%missdata(outcome=cesarean,     label=Caesarean delivery);
%missdata(outcome=ptb,          label=Preterm birth);
%missdata(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 2. MULTIPLE IMPUTATION OF MARITAL STATUS (hypertensive disorders)
  FCS logistic, 50 imputations, seed 20260616. Rubin's rules.
============================================================================*/

data asian_mi;
  set ana.analytic_pooled;
  where asian = 1;

  length bmi5 $7;
  if 13.0 <= bmi <= 69.9 then do;
    if      bmi < 18.5 then bmi5 = 'under';
    else if bmi < 23.0 then bmi5 = 'normal';
    else if bmi < 25.0 then bmi5 = 'hidden';
    else if bmi < 27.5 then bmi5 = 'visover';
    else                    bmi5 = 'obese';
  end;
  else bmi5 = '';
  bmi5_under   = (bmi5 = 'under');
  bmi5_hidden  = (bmi5 = 'hidden');
  bmi5_visover = (bmi5 = 'visover');
  bmi5_obese   = (bmi5 = 'obese');

  length asiangrp $12;
  select (mrace15);
    when (4)  asiangrp = 'AsianIndian';
    when (5)  asiangrp = 'Chinese';
    when (6)  asiangrp = 'Filipino';
    when (7)  asiangrp = 'Japanese';
    when (8)  asiangrp = 'Korean';
    when (9)  asiangrp = 'Vietnamese';
    when (10) asiangrp = 'OtherAsian';
    otherwise asiangrp = '';
  end;

  if dmar not in (1,2) then dmar = .;

  /* complete case on everything except marital status */
  if missing(htn_disorder) or missing(agecat) or missing(meduc)
     or missing(lbo_rec) or missing(precare5) or missing(pay_rec)
     or missing(sex) or missing(smoke) or missing(nativity)
     or missing(year) or bmi5 = '' or asiangrp = '' then delete;
run;

proc mi data=asian_mi nimpute=50 out=asian_imp seed=20260616 noprint;
  class dmar agecat meduc lbo_rec precare5 pay_rec sex smoke nativity
        year asiangrp;
  fcs logistic(dmar = bmi5_under bmi5_hidden bmi5_visover bmi5_obese
                      agecat meduc lbo_rec precare5 pay_rec sex smoke
                      nativity year asiangrp htn_disorder);
  var bmi5_under bmi5_hidden bmi5_visover bmi5_obese
      agecat meduc lbo_rec precare5 pay_rec sex smoke nativity year
      asiangrp htn_disorder dmar;
run;

proc sort data=asian_imp; by _Imputation_; run;

ods rtf exclude all;
proc glimmix data=asian_imp empirical=classical;
  by _Imputation_;
  class dmar(ref='1') agecat(ref='25-29') meduc(ref='3') lbo_rec(ref='1')
        precare5(ref='1') pay_rec(ref='2') sex(ref='F') smoke(ref='0')
        nativity(ref='USborn') year(ref='2016') asiangrp(ref='AsianIndian');
  model htn_disorder = &BMIX asiangrp agecat meduc dmar lbo_rec precare5
        pay_rec sex smoke nativity year
        / dist=poisson link=log solution ddfm=residual;
  random _residual_;
  ods output ParameterEstimates=mi_parms;
run;
ods rtf exclude none;

ods output ParameterEstimates=mi_pooled VarianceInfo=mi_varinfo;
proc mianalyze parms=mi_parms;
  modeleffects &BMIX;
run;

proc sort data=mi_pooled;  by Parm; run;
proc sort data=mi_varinfo; by Parm; run;

data mi_rr;
  merge mi_pooled(keep=Parm Estimate LCLMean UCLMean)
        mi_varinfo(keep=Parm FracMiss);
  by Parm;
  RR    = exp(Estimate);
  RR_lo = exp(LCLMean);
  RR_hi = exp(UCLMean);
run;

title "Section 2 (eTable 6): MI pooled RRs, hypertensive disorders";
proc print data=mi_rr noobs label;
  var Parm RR RR_lo RR_hi FracMiss;
  format RR RR_lo RR_hi 6.2 FracMiss 8.4;
  label RR_lo = 'Lower 95% CL' RR_hi = 'Upper 95% CL'
        FracMiss = 'Fraction missing information';
run;


/*============================================================================
  SECTION 3. HYPERTENSIVE DISORDERS: ECLAMPSIA EXCLUDED; 2022-2023 ONLY
============================================================================*/

title "Section 3a (eTable 6): gestational hypertension only (eclampsia excluded)";
proc glimmix data=asian5 empirical=classical;
  class &CLS asiangrp(ref='AsianIndian');
  model htn_noecl = &BMIX asiangrp &COV / dist=poisson link=log;
  random _residual_;
  %bmi5est;
run;

title "Section 3b (eTable 6): hypertensive disorders, 2022-2023 only";
proc glimmix data=asian5 empirical=classical;
  where year >= 2022;
  class &CLSB dmar(ref='1') year(ref='2022') asiangrp(ref='AsianIndian');
  model htn_disorder = &BMIX asiangrp &COV / dist=poisson link=log;
  random _residual_;
  %bmi5est;
run;


/*============================================================================
  SECTION 4. eTable 5 PANEL A: RACE x BMI SPLINE ROBUSTNESS
  Same restricted cubic spline as 02_main_models.sas, Section 6.
  4a US-born, 4b foreign-born, 4c excluding 2020-2021: 4-df race x BMI test
  4d three-way race x BMI x nativity: 4-df test
============================================================================*/

proc univariate data=wa noprint;
  var bmi;
  output out=_knots pctlpts=5 27.5 50 72.5 95 pctlpre=k;
run;

data _null_;
  set _knots;
  call symputx('k1', k5);
  call symputx('k2', k27_5);
  call symputx('k3', k50);
  call symputx('k4', k72_5);
  call symputx('k5', k95);
run;

data wa_rcs;
  set wa;
  array _k{5} _temporary_ (&k1 &k2 &k3 &k4 &k5);
  array s{4} s1-s4;
  s1 = bmi - 25;
  do _j = 1 to 3;
    s{_j+1} = ( max(bmi-_k{_j},0)**3
              - max(bmi-_k{4},0)**3 * (_k{5}-_k{_j}) / (_k{5}-_k{4})
              + max(bmi-_k{5},0)**3 * (_k{4}-_k{_j}) / (_k{5}-_k{4}) )
              / (_k{5}-_k{1})**2;
  end;
  array a{4}  a_s1-a_s4;     /* race x BMI            */
  array f{4}  f_s1-f_s4;     /* BMI x nativity        */
  array af{4} af_s1-af_s4;   /* race x BMI x nativity */
  do _j = 1 to 4;
    a{_j}  = asian * s{_j};
    f{_j}  = fb * s{_j};
    af{_j} = asian * fb * s{_j};
  end;
  a_fb = asian * fb;
  drop _j;
run;

%macro rcs_test(outcome=, where=1, tag=, label=, strat=0);
  title "Section &tag (eTable 5): race x BMI spline test - &label";
  ods select Contrasts;
  proc glimmix data=wa_rcs empirical=classical;
    where &where;
    %if &strat = 1 %then %do;
      class &CLSNN;
      model &outcome = asian s1-s4 a_s1-a_s4 &COVNN / dist=poisson link=log;
    %end;
    %else %do;
      class &CLS;
      model &outcome = asian s1-s4 a_s1-a_s4 &COV / dist=poisson link=log;
    %end;
    random _residual_;
    contrast "Race x BMI (4 df)" a_s1 1, a_s2 1, a_s3 1, a_s4 1 / chisq;
  run;
  ods select all;
%mend rcs_test;

%macro rcs_3way(outcome=, label=);
  title "Section 4d (eTable 5): race x BMI x nativity (4 df) - &label";
  ods select Contrasts;
  proc glimmix data=wa_rcs empirical=classical;
    class &CLS;
    model &outcome = asian s1-s4 a_s1-a_s4 f_s1-f_s4 a_fb af_s1-af_s4 &COV
          / dist=poisson link=log;
    random _residual_;
    contrast "Race x BMI x nativity (4 df)"
             af_s1 1, af_s2 1, af_s3 1, af_s4 1 / chisq;
  run;
  ods select all;
%mend rcs_3way;

%macro run_spline(outcome=, label=);
  %rcs_test(outcome=&outcome, where=nativity='USborn',       tag=4a, label=&label - US-born, strat=1);
  %rcs_test(outcome=&outcome, where=nativity='Foreignborn',  tag=4b, label=&label - foreign-born, strat=1);
  %rcs_test(outcome=&outcome, where=year not in (2020 2021), tag=4c, label=&label - excluding 2020-2021);
  %rcs_3way(outcome=&outcome, label=&label);
%mend run_spline;

%run_spline(outcome=htn_disorder, label=Hypertensive disorders);
%run_spline(outcome=cesarean,     label=Caesarean delivery);
%run_spline(outcome=ptb,          label=Preterm birth);
%run_spline(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 5. eTable 6 PANEL C: E-VALUES (VanderWeele and Ding, 2017)
  E = RR + sqrt(RR x (RR - 1)), for the point estimate and the lower
  confidence limit. RRs are the Table 3 hypertensive-disorder estimates.
============================================================================*/

data evalues;
  input category $12. rr lcl;
  e_point = rr  + sqrt(rr  * (rr  - 1));
  e_lcl   = lcl + sqrt(lcl * (lcl - 1));
  datalines;
23.0-24.9   1.44 1.40
>=27.5      2.79 2.73
;
run;

title "Section 5 (eTable 6): E-values, hypertensive disorders";
proc print data=evalues noobs;
  format rr lcl e_point e_lcl 6.2;
run;

title;
ods rtf close;
