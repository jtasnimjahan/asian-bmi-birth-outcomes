/******************************************************************************
*** 02_main_models.sas
***
*** Project : Hidden overweight and adverse birth outcomes among Asian mothers
***           in the United States, 2016-2023
*** Purpose : Main models. Modified Poisson regression (log link, empirical
***           sandwich variance) in PROC GLIMMIX.
*** Input   : ana.wa (from 01_data_prep.sas), ana.analytic_pooled
*** Produces: Section 0  Table 2 (Asian vs White crude and adjusted RRs)
***           Section 1  eTable 3 Panel A (complete-case N and events)
***           Section 2  Table 3 and Figure 2 (standardized risks, RDs, RRs)
***           Section 3  LBW RR without subgroup adjustment (Results text)
***           Section 4  Table 4 (subgroup RRs and heterogeneity tests)
***           Section 5  eTable 4 (slope change at BMI 23.0)
***           Section 6  eTable 5 Panel B and overall Panel A row (spline)
***           Section 7  eFigure 2 data (efig2_spline_grid.csv)
***           Section 8  eTable 1 and eFigure 3 values
*** Software: SAS 9.4
*** Author  : Jabin Tasnim Jahan, University at Albany, SUNY
*** Contact : jtasnimjahan@albany.edu
*** Version : 1.0 (November 2026), code for the submitted manuscript
******************************************************************************/

options nocenter ls=200 ps=max mprint;

/* Set these two paths before running */
libname ana "/path/to/your/data";
%let OUTDIR = /path/to/your/output;

ods rtf file="&OUTDIR./02_main_models_results.rtf" style=journal;

%let COV  = agecat meduc dmar lbo_rec precare5 pay_rec sex smoke nativity year;
%let CLS  = agecat(ref='25-29') meduc(ref='3') dmar(ref='1') lbo_rec(ref='1')
            precare5(ref='1') pay_rec(ref='2') sex(ref='F') smoke(ref='0')
            nativity(ref='USborn') year(ref='2016');
%let BMIX = bmi5_under bmi5_hidden bmi5_visover bmi5_obese;

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
  SECTION 0. TABLE 2: ASIAN VS WHITE RISK RATIOS
  0a  n/N (%) by race
  0b  Crude RR (race only)
  0c  Adjusted RR: standard WHO BMI category and all covariates
  0d  Footnote: adjusted RR with missing marital status as a category
============================================================================*/

title "Section 0a (Table 2): events by race";
proc freq data=wa;
  tables raceeth*(htn_disorder cesarean ptb lbw) / nocol nopercent;
run;

%macro rr_disparity(outcome=, label=);
  title "Section 0b (Table 2): crude Asian vs White RR - &label";
  proc glimmix data=wa empirical=classical;
    model &outcome = asian / dist=poisson link=log;
    random _residual_;
    estimate "Asian vs White RR, crude" asian 1 / exp cl;
  run;

  title "Section 0c (Table 2): adjusted Asian vs White RR - &label";
  proc glimmix data=wa empirical=classical;
    class bmi_std(ref='normal') &CLS;
    model &outcome = asian bmi_std &COV / dist=poisson link=log;
    random _residual_;
    estimate "Asian vs White RR, adjusted" asian 1 / exp cl;
  run;

  title "Section 0d (Table 2 footnote): missing marital status as a category - &label";
  proc glimmix data=wa empirical=classical;
    class bmi_std(ref='normal') agecat(ref='25-29') meduc(ref='3') dmar_mi(ref='1')
          lbo_rec(ref='1') precare5(ref='1') pay_rec(ref='2') sex(ref='F')
          smoke(ref='0') nativity(ref='USborn') year(ref='2016');
    model &outcome = asian bmi_std agecat meduc dmar_mi lbo_rec precare5 pay_rec
                     sex smoke nativity year / dist=poisson link=log;
    random _residual_;
    estimate "Asian vs White RR, missing-indicator" asian 1 / exp cl;
  run;
%mend rr_disparity;

%rr_disparity(outcome=htn_disorder, label=Hypertensive disorders);
%rr_disparity(outcome=cesarean,     label=Caesarean delivery);
%rr_disparity(outcome=ptb,          label=Preterm birth);
%rr_disparity(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 1. eTable 3 PANEL A: COMPLETE-CASE N AND EVENTS BY BMI CATEGORY
  Unknown codes (meduc=9, precare5=5, pay_rec=9, lbo_rec=9) are kept as
  categories; only SAS missing values cause complete-case exclusion.
============================================================================*/

%macro nevents(outcome=);
  data _cc;
    set asian5;
    if cmiss(&outcome, agecat, meduc, dmar, lbo_rec, precare5,
             pay_rec, sex, smoke, nativity, year) = 0;
  run;
  title "Section 1: Complete-case N and events by BMI category - &outcome";
  proc freq data=_cc;
    tables bmi5*&outcome / nocol nopercent;
  run;
%mend nevents;

%nevents(outcome=htn_disorder);
%nevents(outcome=cesarean);
%nevents(outcome=ptb);
%nevents(outcome=lbw);


/*============================================================================
  SECTION 2. TABLE 3 AND FIGURE 2
  Subgroup-adjusted model. Marginal standardization over the complete-case
  Asian sample; 95% CIs by the delta method with the sandwich covariance.
============================================================================*/

%macro stdrisk(outcome=, label=);

  proc datasets lib=work nolist;
    delete _cc _des _pe _cov _b _g stdadj_&outcome;
  quit;

  data _cc;
    set asian5;
    if cmiss(&outcome, agecat, meduc, dmar, lbo_rec, precare5,
             pay_rec, sex, smoke, nativity, year, asiangrp) = 0;
  run;

  proc glmselect data=_cc outdesign(addinputvars)=_des noprint;
    class &CLS asiangrp(ref='AsianIndian') / param=ref;
    model &outcome = &BMIX asiangrp &COV / selection=none;
  run;

  %let xlist = &_GLSMOD;
  %let p  = %sysfunc(countw(&xlist, %str( )));
  %let pp = %eval(&p + 1);

  ods rtf exclude all;
  ods output ParameterEstimates=_pe CovB=_cov;
  proc glimmix data=_des empirical=classical;
    model &outcome = &xlist / dist=poisson link=log s covb;
    random _residual_;
  run;
  ods rtf exclude none;

  data _pe; set _pe; where upcase(Effect) not in ('SCALE','RESIDUAL'); run;
  proc transpose data=_pe out=_b prefix=b; var Estimate; run;

  /* Mean predicted risk and gradient under each BMI scenario:
     1 = normal (ref), 2 = under, 3 = hidden, 4 = visible overweight, 5 = obese */
  data _g(keep=scen risk g1-g&pp);
    if _n_ = 1 then set _b(keep=b1-b&pp);
    set _des(keep=&xlist) end=last;
    array x{&p}  &xlist;
    array b{&pp} b1-b&pp;
    array g{&pp} g1-g&pp;
    array S{5,&pp} _temporary_;
    retain _nobs 0;
    if _n_ = 1 then do _s = 1 to 5; do _j = 1 to &pp; S{_s,_j} = 0; end; end;
    _nobs + 1;
    do _s = 1 to 5;
      x{1} = (_s = 2); x{2} = (_s = 3); x{3} = (_s = 4); x{4} = (_s = 5);
      _eta = b{1};
      do _j = 1 to &p; _eta = _eta + b{_j+1} * x{_j}; end;
      _pr = exp(_eta);
      S{_s,1} = S{_s,1} + _pr;
      do _j = 1 to &p; S{_s,_j+1} = S{_s,_j+1} + _pr * x{_j}; end;
    end;
    if last then do _s = 1 to 5;
      scen = _s;
      do _j = 1 to &pp; g{_j} = S{_s,_j} / _nobs; end;
      risk = g{1};
      output;
    end;
  run;

  proc iml;
    use _cov; read all var ("Col1":"Col&pp") into V; close _cov;
    use _g;   read all var ("g1":"g&pp") into G; read all var {risk} into R; close _g;
    lab = {"Normal 18.5-22.9 (ref)", "Underweight <18.5", "Hidden 23.0-24.9",
           "Visible OW 25.0-27.4", "Obese >=27.5"};
    z = 1.959964;
    seR   = sqrt(vecdiag(G * V * G`));
    D     = G[2:5,] - repeat(G[1,], 4, 1);
    RD    = R[2:5] - R[1];
    seRD  = sqrt(vecdiag(D * V * D`));
    L     = G[2:5,] / repeat(R[2:5], 1, &pp) - repeat(G[1,] / R[1], 4, 1);
    lRR   = log(R[2:5] / R[1]);
    seLRR = sqrt(vecdiag(L * V * L`));
    risk_pct = 100*R;  risk_lo = 100*(R - z*seR);  risk_hi = 100*(R + z*seR);
    rd_pp  = {.} // 100*RD;
    rd_lo  = {.} // 100*(RD - z*seRD);
    rd_hi  = {.} // 100*(RD + z*seRD);
    srr    = {1} // exp(lRR);
    srr_lo = {.} // exp(lRR - z*seLRR);
    srr_hi = {.} // exp(lRR + z*seLRR);
    outcome  = repeat("&outcome", 5, 1);
    category = lab`;
    create stdadj_&outcome var {outcome category risk_pct risk_lo risk_hi
                                rd_pp rd_lo rd_hi srr srr_lo srr_hi};
    append; close stdadj_&outcome;
  quit;

  title "Section 2 (Table 3): standardized risk (%), RD (percentage points), RR - &label";
  proc print data=stdadj_&outcome noobs;
    format risk_pct risk_lo risk_hi rd_pp rd_lo rd_hi srr srr_lo srr_hi 8.2;
  run;
%mend stdrisk;

%stdrisk(outcome=htn_disorder, label=Hypertensive disorders);
%stdrisk(outcome=cesarean,     label=Caesarean delivery);
%stdrisk(outcome=ptb,          label=Preterm birth);
%stdrisk(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 3. LOW BIRTHWEIGHT RR WITHOUT SUBGROUP ADJUSTMENT (Results text)
  Same complete-case sample and covariates as Table 3, asiangrp omitted.
============================================================================*/

title "Section 3: Low birthweight BMI RRs without subgroup adjustment";
proc glimmix data=asian5 empirical=classical;
  class &CLS;
  model lbw = &BMIX &COV / dist=poisson link=log;
  random _residual_;
  %bmi5est;
run;


/*============================================================================
  SECTION 4. TABLE 4: BMI 23.0-24.9 RR BY DETAILED ASIAN SUBGROUP
  4a  Five-category BMI model fitted within each subgroup
  4b  6-df Wald test of BMI category x subgroup, mothers with BMI 18.5-24.9
============================================================================*/

%macro subgroups;
  %local g o grp out;
  %let grp = AsianIndian Chinese Filipino Japanese Korean Vietnamese OtherAsian;
  %let out = htn_disorder cesarean ptb lbw;
  %do g = 1 %to 7;
    %do o = 1 %to 4;
      title "Section 4a (Table 4): %scan(&grp,&g) - %scan(&out,&o)";
      proc glimmix data=asian5 empirical=classical;
        where asiangrp = "%scan(&grp,&g)";
        class &CLS;
        model %scan(&out,&o) = &BMIX &COV / dist=poisson link=log;
        random _residual_;
        %bmi5est;
      run;
    %end;
  %end;
%mend subgroups;
%subgroups;

data asian_nh;
  set asian5;
  where bmi5 in ('normal','hidden');
run;

%macro hetero(outcome=, label=);
  title "Section 4b (Table 4): heterogeneity across subgroups, 6 df - &label";
  proc glimmix data=asian_nh empirical=classical;
    class bmi5(ref='normal') asiangrp(ref='AsianIndian') &CLS;
    model &outcome = bmi5 asiangrp bmi5*asiangrp &COV
          / dist=poisson link=log chisq;
    random _residual_;
  run;
%mend hetero;

%hetero(outcome=htn_disorder, label=Hypertensive disorders);
%hetero(outcome=cesarean,     label=Caesarean delivery);
%hetero(outcome=ptb,          label=Preterm birth);
%hetero(outcome=lbw,          label=Low birthweight);
/* Report the bmi5*asiangrp row of the Type III table. */


/*============================================================================
  SECTION 5. eTable 4: SLOPE CHANGE AT BMI 23.0
  Asian mothers with BMI 18.5-27.4. exp(hinge23) = ratio of the RR per
  1 kg/m2 above versus below 23.0.
============================================================================*/

%macro hinge(outcome=, label=);
  title "Section 5 (eTable 4): slope change at BMI 23.0 - &label";
  proc glimmix data=asian5 empirical=classical;
    where 18.5 <= bmi < 27.5;
    class &CLS asiangrp(ref='AsianIndian');
    model &outcome = bmi hinge23 asiangrp &COV / dist=poisson link=log;
    random _residual_;
    estimate "RR per 1 kg/m2, BMI <23"   bmi 1           / exp cl;
    estimate "RR per 1 kg/m2, BMI >=23"  bmi 1 hinge23 1 / exp cl;
    estimate "Ratio of per-unit RRs"     hinge23 1       / exp cl;
  run;
%mend hinge;

%hinge(outcome=htn_disorder, label=Hypertensive disorders);
%hinge(outcome=cesarean,     label=Caesarean delivery);
%hinge(outcome=ptb,          label=Preterm birth);
%hinge(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 6. eTable 5: RACE x BMI WITH A RESTRICTED CUBIC SPLINE
  Harrell RCS, 5 knots at the 5th, 27.5th, 50th, 72.5th and 95th percentiles
  of BMI in the Asian + White sample. s1 = BMI - 25 (linear), s2-s4
  nonlinear. a_s1-a_s4 = asian x s1-s4. 4-df test of a_s1-a_s4 and
  Asian vs White RRs at the same BMI.
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
%put NOTE: RCS knots = &k1 &k2 &k3 &k4 &k5;

/* Basis macro: used once per data step (it declares arrays) */
%macro rcs(x);
  _s1 = &x - 25;
  array _k{5} _temporary_ (&k1 &k2 &k3 &k4 &k5);
  array _s{2:4} _s2-_s4;
  _den = (_k{5} - _k{4});
  do _j = 1 to 3;
    _s{_j+1} = ( max(&x-_k{_j},0)**3
               - max(&x-_k{4},0)**3 * (_k{5}-_k{_j}) / _den
               + max(&x-_k{5},0)**3 * (_k{4}-_k{_j}) / _den )
               / (_k{5}-_k{1})**2;
  end;
%mend rcs;

data wa_rcs;
  set wa;
  %rcs(bmi);
  s1 = _s1; s2 = _s2; s3 = _s3; s4 = _s4;
  a_s1 = asian*s1; a_s2 = asian*s2; a_s3 = asian*s3; a_s4 = asian*s4;
  drop _s1-_s4 _den _j;
run;

%let BMIGRID = 20 23 25 27.5 30 35;
data _null_;
  do _i = 1 to %sysfunc(countw(&BMIGRID, %str( )));
    b = input(scan("&BMIGRID", _i, ' '), best.);
    %rcs(b);
    call symputx(cats('c', _i, '_1'), _s1);
    call symputx(cats('c', _i, '_2'), _s2);
    call symputx(cats('c', _i, '_3'), _s3);
    call symputx(cats('c', _i, '_4'), _s4);
  end;
run;

%macro race_rcs(outcome=, label=);
  %local i v;
  title "Section 6 (eTable 5): race x BMI, restricted cubic spline - &label";
  proc glimmix data=wa_rcs empirical=classical;
    class &CLS;
    model &outcome = asian s1-s4 a_s1-a_s4 &COV / dist=poisson link=log s;
    random _residual_;
    contrast "Race x BMI interaction (4 df)"
             a_s1 1, a_s2 1, a_s3 1, a_s4 1 / chisq;
    %do i = 1 %to %sysfunc(countw(&BMIGRID, %str( )));
      %let v = %scan(&BMIGRID, &i, %str( ));
      estimate "Asian vs White, BMI &v"
               asian 1 a_s1 &&c&i._1 a_s2 &&c&i._2 a_s3 &&c&i._3 a_s4 &&c&i._4
               / exp cl;
    %end;
  run;
%mend race_rcs;

%race_rcs(outcome=htn_disorder, label=Hypertensive disorders);
%race_rcs(outcome=cesarean,     label=Caesarean delivery);
%race_rcs(outcome=ptb,          label=Preterm birth);
%race_rcs(outcome=lbw,          label=Low birthweight);


/*============================================================================
  SECTION 7. eFigure 2 DATA
  RR at each BMI from 17 to 40 (step 0.5), relative to White mothers at
  BMI 25. Same model as Section 6. Output: efig2_spline_grid.csv
============================================================================*/

%let GLO = 17;
%let GHI = 40;
%let GBY = 0.5;

data _null_;
  %rcs(25);
  call symputx('r1', _s1);
  call symputx('r2', _s2);
  call symputx('r3', _s3);
  call symputx('r4', _s4);
run;

filename estf temp;
data _null_;
  file estf lrecl=2000;
  length lab $12;
  r1 = &r1; r2 = &r2; r3 = &r3; r4 = &r4;
  do b = &GLO to &GHI by &GBY;
    %rcs(b);
    d1 = _s1 - r1; d2 = _s2 - r2; d3 = _s3 - r3; d4 = _s4 - r4;
    if abs(b - 25) > 1e-8 then do;
      lab = cats('W_', put(b, 4.1));
      put 'estimate "' lab +(-1) '" s1 ' d1 best16. ' s2 ' d2 best16.
          ' s3 ' d3 best16. ' s4 ' d4 best16. ' / exp cl;';
    end;
    lab = cats('A_', put(b, 4.1));
    put 'estimate "' lab +(-1) '" asian 1 s1 ' d1 best16. ' s2 ' d2 best16.
        ' s3 ' d3 best16. ' s4 ' d4 best16.
        ' a_s1 ' _s1 best16. ' a_s2 ' _s2 best16.
        ' a_s3 ' _s3 best16. ' a_s4 ' _s4 best16. ' / exp cl;';
  end;
run;

%macro grid(outcome=);
  proc datasets lib=work nolist; delete _est; quit;
  ods select Estimates;
  ods output Estimates=_est;
  proc glimmix data=wa_rcs empirical=classical;
    class &CLS;
    model &outcome = asian s1-s4 a_s1-a_s4 &COV / dist=poisson link=log;
    random _residual_;
    %include estf;
  run;
  ods select all;

  data _est_&outcome(keep=outcome group bmi rr lo hi);
    length outcome $12 group $5;
    set _est;
    outcome = "&outcome";
    group   = ifc(substr(Label,1,1) = 'A', 'Asian', 'White');
    bmi     = input(substr(Label,3), best.);
    rr = ExpEstimate; lo = ExpLower; hi = ExpUpper;
  run;
%mend grid;

%grid(outcome=htn_disorder);
%grid(outcome=cesarean);
%grid(outcome=ptb);
%grid(outcome=lbw);

data efig2;
  set _est_htn_disorder _est_cesarean _est_ptb _est_lbw;
run;

proc export data=efig2 outfile="&OUTDIR./efig2_spline_grid.csv"
            dbms=csv replace;
run;


/*============================================================================
  SECTION 8. eTable 1 AND eFigure 3
  8a  eFigure 3: % with BMI 23.0-24.9 among mothers with BMI 18.5-24.9,
      by age group, all Asian and Asian Indian (mrace15 = 4)
  8b  eTable 1: Asian Indian mothers, outcome-specific N and events by
      BMI group (hidden vs normal) and age group
============================================================================*/

data std_normal;
  set ana.analytic_pooled(keep=asian mrace15 bmi agecat);
  where asian = 1 and 18.5 <= bmi < 25.0;
  length group $12;
  hidden_pct = 100 * (bmi >= 23.0);
  group = 'All Asian';    output;
  if mrace15 = 4 then do; group = 'AsianIndian'; output; end;
run;

title "Section 8a (eFigure 3): % with BMI 23.0-24.9 among standard-normal mothers";
proc means data=std_normal mean n maxdec=1 nolabels;
  class group agecat;
  types group group*agecat;
  var hidden_pct;
run;

data ai;
  set ana.analytic_pooled(keep=asian mrace15 bmi agecat htn_disorder cesarean ptb lbw);
  where asian = 1 and mrace15 = 4 and 18.5 <= bmi < 25.0;
  length bmigrp $6;
  if bmi >= 23.0 then bmigrp = 'hidden'; else bmigrp = 'normal';
run;

title "Section 8b (eTable 1): outcome-specific N and events, Asian Indian mothers";
proc means data=ai n sum mean maxdec=4 nolabels;
  class bmigrp agecat;
  types bmigrp bmigrp*agecat;
  var htn_disorder cesarean ptb lbw;
run;

title;
ods rtf close;
