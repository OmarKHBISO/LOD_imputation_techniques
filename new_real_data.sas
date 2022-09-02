libname realdata '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS';
libname sasdata '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS';
%let out124=\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\;



/* Compilation of macro "Macro_AuC2" */
%include '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\Clinique\PHN3314\18.0124.P\Programmes\3314_180124_Macro_AuC2.sas';


%let LoD = 0.7;

/* ================================================================ */
/* 			1- Data 												*/
/* ================================================================ */

proc import datafile="\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\Clinique\PN3323_DAP Canine Line CN\2020027_CPV_CDV_CAV1_CAV2\Datasets\2020027-shedding-challenge-post-V4-for stats-postAQ_QCed.xlsx" 
	out=V4_value_wide replace dbms=xlsx;
	range="Values$A1:L21";
	getnames=yes;
run;
proc import datafile="\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\Clinique\PN3323_DAP Canine Line CN\2020027_CPV_CDV_CAV1_CAV2\Datasets\2020027-shedding-challenge-post-V4-for stats-postAQ_QCed.xlsx" 
	out=V4_censored_wide replace dbms=xlsx;
	range="Censored$A1:L21";
	getnames=yes;
run;
proc import datafile="\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\Clinique\PN3323_DAP Canine Line CN\2020027_CPV_CDV_CAV1_CAV2\Datasets\Transposed-2020027-shedding-challenge-for stats.xlsx" 
	out=V3_value replace dbms=xlsx;
	sheet="Value-postV3-transposed";
	getnames=yes;
run;
proc import datafile="\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\Clinique\PN3323_DAP Canine Line CN\2020027_CPV_CDV_CAV1_CAV2\Datasets\Transposed-2020027-shedding-challenge-for stats.xlsx" 
	out=V3_censored replace dbms=xlsx;
	sheet="Censored-postV3-transposed";
	getnames=yes;
run;

/*=========================== data post V3 =============================*/
proc sql;
create table shedding_V3 as
  select a.Group, a.Animal, a.Chrono, a.Censored, b.Result
  from V3_censored as a
  left join V3_value as b 
    on a.Group = b.Group 
    and a.Animal = b.Animal
    and a.Chrono = b.Chrono
;
quit;


/*================= impute velues below LOD by LOD/2 and (1/sqrt 2)*LOD ============== */

data shedding_V3;
set shedding_V3 (keep=group animal chrono result censored);
Time= input (compress (CHrono,"T"), best2.);
Group_ID=left(trim(Group)) || "-" || left(trim(Animal));
if censored in (1,2) and result=0.7 then titer_imput= (Result/2); else titer_imput=Result;
if censored in (1,2) and result=0.7 then titer_imput2=((1/sqrt(2))*(Result)); else titer_imput2=Result;
run;

/*=========================== data post V3 ================================*/
proc transpose data=V4_value_wide out=V4_value_long;
	by Group Animal;
run;
data V4_value;
	set V4_value_long;
	rename _label_=Chrono COL1=Result;
run;
proc transpose data=V4_censored_wide out=V4_censored_long;
	by Group Animal;
run;
data V4_censored;
	set V4_censored_long;
	rename _label_=Chrono COL1=Censored;
run;
proc sql;
create table shedding_V4 as
  select a.Group, a.Animal, a.Chrono, a.Censored, b.Result
  from V4_censored as a
  left join V4_value as b 
    on a.Group = b.Group 
    and a.Animal = b.Animal
    and a.Chrono = b.Chrono
;
quit;
/*================= impute velues below LOD by LOD/2 and (1/sqrt 2)*LOD ============== */
data shedding_V4;
set shedding_V4 (keep=group animal chrono result censored);
Time= input (compress (CHrono,"T"), best2.);
Group_ID=left(trim(Group)) || "-" || left(trim(Animal));
if censored in (1,2) and result=0.7 then titer_imput= Result/2; else titer_imput=Result;
if censored in (1,2) and result=0.7 then titer_imput2=((1/sqrt(2))*(Result)); else titer_imput2=Result;
run;


/*===================================== titer ==============================*/
proc sort data=shedding_V3; by group animal time;run;
proc sort data=shedding_v4; by group animal time; run;

/*=================== so good sgpannel plot =====================*/

proc sgpanel data=shedding_v3;*where (Group="A") ;
	panelby Group_id/ novarname rows=4 columns=5 ;
	refline 0.7 / label='LoD' labelpos=min lineattrs=(color=green);
	series x=Time y=Result / group=Group transparency=0 markers lineattrs=(pattern=dot);
	keylegend / title="Animal ID";
	rowaxis label="Titer";
	colaxis label="Time (days)";
run;
proc print data=shedding_V3;run;
/*============================ AUC Computation ===============================*/
proc sort data=shedding_V3; by Group Animal Time; run;


/*======================= PAIRWISE COMPARISON ==================================*/

data shedding_V3_new; set shedding_V3; where (group ="A") OR (group= "D");

if time= -1 then do; t1=1 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 2 then do; t1=0 ; t2=1 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 3 then do; t1=0 ; t2=0 ; t3=1 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 4 then do; t1=0 ; t2=0 ; t3=0 ; t4=1 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 5 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=1 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 6 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=1 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 7 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=1; t8=0 ; t9=0 ; t10=0;
end;
 if time= 8 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=1 ; t9=0 ; t10=0;
end;
 if time= 9 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=1 ; t10=0;
end;
 if time= 10 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=1;
end;
if Group= "A"  then do;  trtA=1; Control=0; end;
if Group= "D" then do;  trtA=0; Control=1; end;
/*interaction trtA * time*****/
int1= trtA*t1 ; int2= trtA*t2 ; int3= trtA*t3; int4= trtA*t4; int5= trtA*t5; int6= trtA*t6; int7= trtA*t7; int8=trtA*t8; int9=trtA*t9 ; int10=trtA*t10;
/* interaction trt B * time **********/
int11= Control*t1 ; int12= Control*t2 ; int13= Control*t3; int14= Control*t4; int15= Control*t5; int16= Control*t6; int17= Control*t7; int18=Control*t8; int19=Control*t9 ; int20=Control*t10;
run;

proc print data=shedding_V3_new;run;

/*====================================================================================*/
/****************Apply the Nlmixed procedure  Comparison of A versus C*****************/
/*====================================================================================*/
ods output additionalestimates= Pdiff_NlA;
ods output fitstatistics= fit_NlA;
/*ods output ConvergenceStatus=Convergence_NL&c;*/
proc nlmixed data=shedding_V3_new; 

bounds sigsq1 sigsqe >0;
pi= 2*arsin(1);
mu= beta1*trtA + beta2*Control + beta3*t1 + beta4*t2 + beta5*t3 + beta6*t4 + beta7*t5 + beta8*t6 + beta9*t7 + beta10*t8 + beta11*t9 + beta12*t10 +
beta13*int1 + beta14*int2 + beta15*int3 + beta16*int4 + beta17*int5 + beta18*int6 +beta19*int7 + beta20*int8 + beta21*int9 + beta22*int10 + beta23*int11 + beta24*int12
				+beta25*int13+beta26*int14+beta27*int15+beta28*int16+ beta29*int17+ beta30*int18+ beta31*int19+ beta32*int20+ai;
if titer_imput ne . then ll = (1/sqrt(2*pi*sigsqe))*exp(-(titer_imput-mu)**2/(2*sigsqe));
if titer_imput = . then ll = probnorm((titer_imput-mu)/sqrt(sigsqe));
	L = log(ll);
model titer_imput ~ general (L);

random ai ~ normal(0,sigsq1) subject=Group_id;
estimate 'vaccineA_control_diff' beta1*(9) + beta2*(-9) + beta13*(0.5)+ beta14*(1)+ beta15*(1) + beta16*(1) + beta17*(1) + beta18*(1) + beta19*(1) + beta20*(1) + beta21*(1) + beta22*(0.5)+ beta23*(-0.5) + beta24*(-1) + beta25*(-1)+beta26*(-1)+ beta27*(-1)+beta28*(-1)+beta29*(-1)+beta30*(-1)+beta31*(-1)+beta32*(-0.5);

run; 
/*==============================================================================================*/
/*******************************Apply the Nlmixed procedure for B versus C *********************/
/*=============================================================================================*/
/*======================= PAIRWISE COMPARISON B versus D ==================================*/

data shedding_V3_new; set shedding_V3; where (group ="B") OR (group= "D");

if time= -1 then do; t1=1 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 2 then do; t1=0 ; t2=1 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 3 then do; t1=0 ; t2=0 ; t3=1 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 4 then do; t1=0 ; t2=0 ; t3=0 ; t4=1 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 5 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=1 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 6 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=1 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 7 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=1; t8=0 ; t9=0 ; t10=0;
end;
 if time= 8 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=1 ; t9=0 ; t10=0;
end;
 if time= 9 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=1 ; t10=0;
end;
 if time= 10 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=1;
end;
if Group= "B"  then do;  trt=1; Control=0; end;
if Group= "D" then do;  trt=0; Control=1; end;
/*interaction trt * time*****/
int1= trt*t1 ; int2= trt*t2 ; int3= trt*t3; int4= trt*t4; int5= trt*t5; int6= trt*t6; int7= trt*t7; int8=trt*t8; int9=trt*t9 ; int10=trt*t10;
/* interaction trt B * time **********/
int11= Control*t1 ; int12= Control*t2 ; int13= Control*t3; int14= Control*t4; int15= Control*t5; int16= Control*t6; int17= Control*t7; int18=Control*t8; int19=Control*t9 ; int20=Control*t10;
run;

proc print data=shedding_V3_new;run;

/*====================================================================================*/
/****************Apply the Nlmixed procedure  Comparison of A versus C*****************/
/*====================================================================================*/
ods output additionalestimates= Pdiff_NlB;
ods output Fitstatistics = Fit_NlB;
proc nlmixed data=shedding_V3_new; 

bounds sigsq1 sigsqe >0;
pi= 2*arsin(1);
mu= beta1*trt + beta2*Control + beta3*t1 + beta4*t2 + beta5*t3 + beta6*t4 + beta7*t5 + beta8*t6 + beta9*t7 + beta10*t8 + beta11*t9 + beta12*t10 +
beta13*int1 + beta14*int2 + beta15*int3 + beta16*int4 + beta17*int5 + beta18*int6 +beta19*int7 + beta20*int8 + beta21*int9 + beta22*int10 + beta23*int11 + beta24*int12
				+beta25*int13+beta26*int14+beta27*int15+beta28*int16+ beta29*int17+ beta30*int18+ beta31*int19+ beta32*int20+ai;
if titer_imput ne . then ll = (1/sqrt(2*pi*sigsqe))*exp(-(titer_imput-mu)**2/(2*sigsqe));
if titer_imput = . then ll = probnorm((titer_imput-mu)/sqrt(sigsqe));
	L = log(ll);
model titer_imput ~ general (L);

random ai ~ normal(0,sigsq1) subject=Group_id;
estimate 'vaccineB_control_diff' beta1*(9) + beta2*(-9) + beta13*(0.5)+ beta14*(1)+ beta15*(1) + beta16*(1) + beta17*(1) + beta18*(1) + beta19*(1) + beta20*(1) + beta21*(1) + beta22*(0.5)+ beta23*(-0.5) + beta24*(-1) + beta25*(-1)+beta26*(-1)+ beta27*(-1)+beta28*(-1)+beta29*(-1)+beta30*(-1)+beta31*(-1)+beta32*(-0.5);

run; 








/*==============================================================================================*/
/*******************************PAIRWISE COMPARISON C versus D *******************************************/
/*=============================================================================================*/
/*===============================================================================================*/

data shedding_V3_new; set shedding_V3; where (group ="C") OR (group= "D");

if time= -1 then do; t1=1 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 2 then do; t1=0 ; t2=1 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 3 then do; t1=0 ; t2=0 ; t3=1 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 4 then do; t1=0 ; t2=0 ; t3=0 ; t4=1 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 5 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=1 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 6 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=1 ; t7=0; t8=0 ; t9=0 ; t10=0;
end;
 if time= 7 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=1; t8=0 ; t9=0 ; t10=0;
end;
 if time= 8 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=1 ; t9=0 ; t10=0;
end;
 if time= 9 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=1 ; t10=0;
end;
 if time= 10 then do; t1=0 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=1;
end;
if Group= "D"  then do;  trt=1; Control=0; end;
if Group= "C" then do;  trt=0; Control=1; end;
/*interaction trt * time*****/
int1= trt*t1 ; int2= trt*t2 ; int3= trt*t3; int4= trt*t4; int5= trt*t5; int6= trt*t6; int7= trt*t7; int8=trt*t8; int9=trt*t9 ; int10=trt*t10;
/* interaction trt B * time **********/
int11= Control*t1 ; int12= Control*t2 ; int13= Control*t3; int14= Control*t4; int15= Control*t5; int16= Control*t6; int17= Control*t7; int18=Control*t8; int19=Control*t9 ; int20=Control*t10;
run;

proc print data=shedding_V3_new;run;

/*====================================================================================*/
/****************Apply the Nlmixed procedure  Comparison of A versus C*****************/
/*====================================================================================*/
ods output additionalestimates= Pdiff_NLD;
ods output Fitstatistics=Fit_NLC;
proc nlmixed data=shedding_V3_new; 

bounds sigsq1 sigsqe >0;
pi= 2*arsin(1);
mu= beta1*control + beta2*trt + beta3*t1 + beta4*t2 + beta5*t3 + beta6*t4 + beta7*t5 + beta8*t6 + beta9*t7 + beta10*t8 + beta11*t9 + beta12*t10 +
beta13*int1 + beta14*int2 + beta15*int3 + beta16*int4 + beta17*int5 + beta18*int6 +beta19*int7 + beta20*int8 + beta21*int9 + beta22*int10 + beta23*int11 + beta24*int12
				+beta25*int13+beta26*int14+beta27*int15+beta28*int16+ beta29*int17+ beta30*int18+ beta31*int19+ beta32*int20+ai;
if titer_imput ne . then ll = (1/sqrt(2*pi*sigsqe))*exp(-(titer_imput-mu)**2/(2*sigsqe));
if titer_imput = . then ll = probnorm((titer_imput-mu)/sqrt(sigsqe));
	L = log(ll);
model titer_imput ~ general (L);

random ai ~ normal(0,sigsq1) subject=Group_id;
estimate 'vaccineD_control_diff' beta1*(9) + beta2*(-9) + beta13*(0.5)+ beta14*(1)+ beta15*(1) + beta16*(1) + beta17*(1) + beta18*(1) + beta19*(1) + beta20*(1) + beta21*(1) + beta22*(0.5)+ beta23*(-0.5) + beta24*(-1) + beta25*(-1)+beta26*(-1)+ beta27*(-1)+beta28*(-1)+beta29*(-1)+beta30*(-1)+beta31*(-1)+beta32*(-0.5);

run; 

/*============= using 1/sqrt2 * LOD imputation ====================*/


proc print data=pdiff_nlA;run;

proc print data=pdiff_nlB;run;

proc print data=pdiff_nlC;run;


/*============= exporting results ===================*/

data pdiff_nlA; set pdiff_nlA; method= "LOD/2*Nlmixed";run;
data pdiff_nlB; set pdiff_nlB; method="LOD/2*Nlmixed";run;
data pdiff_nlD; set pdiff_nlD; method="LOD/2*Nlmixed";run;

data sasdata.merging_pdiff; merge pdiff_nlA pdiff_nlB pdiff_nlD; by label;run;

data fit_nlA ; set fit_nlA; Treatment="A versus C"; drop A_versus_C;  run;
data fit_nlB ; set fit_nlB; Treatment="B versus C";drop B_versus_C run;
data fit_nlD ; set fit_nlC; Treatment="D versus C";run;

data fit_nl; merge fit_nlA fit_nlB fit_nlD; by Treatment;run;



ods rtf file='\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\new_merging.doc' style=statistical;
options nodate pageno=1 linesize=80 pagesize=60;	
proc print data=sasdata.merging_pdiff ;
proc print data=fit_nl; 

   title 'Energy Expenditures for Each Region';
   title2 '(millions of dollars)';
run; ods rtf close; 


/*======================= plot differences in mean AUC by group ==========================*/
data sasdata.final_merging ; set sasdata.merging_pdiff;
Estimate2= (-1)*Estimate ;
Lower2= (-1)*Lower;
Upper2= (-1)*Upper;
P_value= Probt;
CI='[ ' || strip(put(round((-1)*Upper,0.01),9.2)) || ' ; ' || strip(put(round((-1)*Lower,0.01),9.2)) || ' ]'; 
drop  Lower Upper DF Probt tValue Estimate Alpha;
run;

/*================== PLOT NLMIXED PROCEDURE =========================*/

proc sgplot data=sasdata.final_merging;
	refline 0 / axis=x lineattrs=(pattern=MediumDash);
	scatter x=Estimate2 y=label / xerrorlower=Lower2 xerrorupper=Upper2;
	xaxis label='AUC Difference Estimate using hybrid imputation technique LOD/2-Nlmixed' values=(0 to 30);
	yaxis label='Comparison between vaccine groups versus control ';
run;

/*======================= Simple imputation technique ============================*/
/*================================================================================*/

proc print data=shedding_V3;run;



/* ================================================================ */
/* 			3.1- AUC computation									*/
/* ================================================================ */

/************** run the macro AUC *******************/
%macro AUC(multdata);
	%let c=1;
	%let dep = %scan(&multdata, &c, ' ');
	%do %until(&dep = end);

		data AUCRep&c;
			set shedding_V3;
			if (time=1|time=10) then do;
				if &dep=. then &dep=0;
			end;
			if (&dep=.) then delete;
			drop LagTime LagValue;
			LagTime = LAG(time);
			LagValue = LAG(&dep);

			if time = 25 then do;
				LagTime = 0;
				LagValue = 0;
				SumTrapezoid&c=0; **My addition: Does AUC for 'by id';
				Trapezoid&c = .;*My addition: Make sure that trapezoid = 0 at first time point;
				SumTrapezoid&c + Trapezoid&c;
			end;
			else do;
				Trapezoid&c = (time-LagTime)*(&dep+LagValue)/2; 
				SumTrapezoid&c + Trapezoid&c;
			end;
		run;
	%let c = %eval(&c + 1);
    %let dep = %scan(&multdata, &c);
	%end;

	proc sort data=AUCRep1;by group_id time;run;
	proc sort data=AUCRep2;by group_id time;run;
	proc sort data=AUCRep3;by group_id time;run;
	proc sort data=AUCRep4;by group_id time;run;

	data AUCRep;
		merge AUCRep1 AUCRep2 AUCRep3 AUCRep3 AUCRep4;
		by group_id time tt;
	run;
	/*proc print data=AUCRep;run;*/
	data AUCBySubjectRep;
		set AUCRep;
		if time = 29;
		AUC1=SumTrapezoid1;AUC2=SumTrapezoid2;AUC3=SumTrapezoid3;AUC4=SumTrapezoid4;
		keep group_id AUC1 AUC2 AUC3 AUC4  Groupc;
	run;

	proc sort data=AUCBySubjectRep; by group_id;run;
	proc sort data=AUCRep;by group_id;run;
	data replacement;
		merge AUCRep AUCBySubjectRep;
		by group_id;
	run;

%mend;
%AUC(multdata=Y1 Y2 Y3 Y4  end); 


