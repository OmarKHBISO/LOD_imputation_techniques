
 libname sasdata '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS';
libname safe '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\LOD_Methods_code\new_datasets';
%GLOBAL n iter multfactor factor ddfm letmesee rho lod;
%PUT _USER_;

options mprint;
%macro NL(NLdata);
	%let d=1;
	%let file = %scan(&NLdata, &d, ' ');
	%do %until(&file = end);






		data nl;
		set safe.&file;where (iter=40);
/*			if conc = . then Ystar1 = (1*(LOD/10)); else Ystar1 = conc;*/
			if conc = . then Ystar2 = ((1/2)*(LOD/10)); else Ystar2 = conc;
			if conc = . then Ystar3 = ((1/sqrt(2))*(LOD/10)); else Ystar3 = conc;
/*			if conc = . then Ystar4 = ((0)*(LOD/10)); else Ystar4 = conc;*/

if time= 1 then do; t1=1 ; t2=0 ; t3=0 ; t4=0 ; t5=0 ; t6=0 ; t7=0; t8=0 ; t9=0 ; t10=0;
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
if trt= "A" then do;  trtA=1; trtB=0; end;
if trt= "B" then do;  trtA=0; trtB=1; end;
/*interaction trtA * time*****/
int1= trtA*t1 ; int2= trtA*t2 ; int3= trtA*t3; int4= trtA*t4; int5= trtA*t5; int6= trtA*t6; int7= trtA*t7; int8=trtA*t8; int9=trtA*t9 ; int10=trtA*t10;
/* interaction trt B * time **********/
int11= trtB*t1 ; int12= trtB*t2 ; int13= trtB*t3; int14= trtB*t4; int15= trtB*t5; int16= trtB*t6; int17= trtB*t7; int18=trtB*t8; int19=trtB*t9 ; int20=trtB*t10;
run;

/*======================================= tenter de nouvelles méthodes ===================================*/











proc sort data=NL;by iter loop; run; 
	data _null_;
			set NL;
			call symput('gamma',left(gamma));
			call symput('lod',left(lod));
			call symput('corr',left(corr));
			call symput('iter',left(iter));
		run;

		
%macro nl_ttest(multdata);
			%let c=1;
			%let dep = %scan(&multdata, &c, ' ');
			%do %until(&dep = end);
				%put In NLMixed: File=&File data=&dep;
		%put _USER_;

		
dm "out;clear;";
ods output additionalestimates= Pdiff_NL&c;
/*ods output ConvergenceStatus=Convergence_NL&c;*/
proc nlmixed data=NL; 

bounds sigsq1 sigsqe >0;
pi= 2*arsin(1);
mu= beta1*trtA + beta2*trtB + beta3*t1 + beta4*t2 + beta5*t3 + beta6*t4 + beta7*t5 + beta8*t6 + beta9*t7 + beta10*t8 + beta11*t9 + beta12*t10 +
beta13*int1 + beta14*int2 + beta15*int3 + beta16*int4 + beta17*int5 + beta18*int6 +beta19*int7 + beta20*int8 + beta21*int9 + beta22*int10 + beta23*int11 + beta24*int12
				+beta25*int13+beta26*int14+beta27*int15+beta28*int16+ beta29*int17+ beta30*int18+ beta31*int19+ beta32*int20+ai;
if &dep ne . then ll = (1/sqrt(2*pi*sigsqe))*exp(-(&dep-mu)**2/(2*sigsqe));
if &dep = . then ll = probnorm((&dep-mu)/sqrt(sigsqe));
	L = log(ll);
model &dep ~ general (L);

random ai ~ normal(0,sigsq1) subject=ID;
estimate 'trt_diff' beta1*(9) + beta2*(-9) + beta13*(0.5)+ beta14*(1)+ beta15*(1) + beta16*(1) + beta17*(1) + beta18*(1) + beta19*(1) + beta20*(1) + beta21*(1) + beta22*(0.5)+ beta23*(-0.5) + beta24*(-1) + beta25*(-1)+beta26*(-1)+ beta27*(-1)+beta28*(-1)+beta29*(-1)+beta30*(-1)+beta31*(-1)+beta32*(-0.5);
by loop iter;
run;

/*data safe.ConvergenceNL&c;*/
/*			set Convergence_NL&c;*/
/*			format _all_;*/
/*		run;*/



                ods output "Statistics"=lsmeans_ttest5;
				ods output "T-Tests"=pdiff_ttest5;
				proc ttest data=nl;
					class trt;
					var AUC;
					by loop iter;
					where  time = 1;
				run;


				data lsmeans_ttest25 ;
					set lsmeans_ttest5;
					where Class= 'Diff (1-2)';
					keep loop iter LowerCLMean UpperCLMean mean StdErr;
					format _all_;
				run;
dm "out;clear;log;clear;";

				data pdiff_ttest25 ;
					set pdiff_ttest5;
					where Method='Satterthwaite';
					keep loop iter df Probt;
					format _all_;
				run;
dm "out;clear;log;clear;";
		
				data test_ttest5;
					merge pdiff_ttest25 lsmeans_ttest25;
					by loop ;
					iter=iter;
					method="cc";
					if ProbT <= 0.05 then rejectNull_1 = 1;
					else rejectNull_1 = 0;
					if ProbT = . then rejectNull_1 = .; 	format _all_;
run;

                dm "out;clear;log;clear;";
				proc sort data=test_ttest5; by iter loop;run;
/******************things are done here ***************************/



	proc means data=test_ttest5 maxdec=2;where (iter=&iter);
	var rejectNull_1 Mean  Stderr ;
	output out=mean_mixed5 mean=power_mixed5  Mean5 Stderr5 ;
	by iter;
	run; 

	data test_NL&c;
			set pdiff_NL&c;
			iter=&iter;
			lod=&lod;
			method="nl";
			if ProbT <= 0.05 then rejectNull = 1;
			else rejectNull = 0;
			if ProbT = . then rejectNull = .;
			pNL=Probt;
			Stderr = StandardError;
			 		run;

proc sort data= test_ttest5 out=test_ttest5 nodupkey; by loop; run;

data nlmixed_bias&&c; merge test_nl&c test_ttest5; by loop iter ;run;


data nlmixed_bias&c; set nlmixed_bias&c;
	Bias_NL = Estimate - mean;
	method_nl =&c;

	
run;





	






		proc print data=nlmixed_bias&c;run;

		proc means data=nlmixed_bias&c	 maxdec=2;
		var rejectNull_1 rejectnull Estimate mean Stderr Bias_NL lod method_nl;
	output out=mean_NL&c mean= power_cc power_nl&c  meanEstimate_NL meanEstimate_cc meanStderr_NL meanBias_NL lod method_nl;
		by iter;
run;

	data power_nlmixed&c; set mean_NL&c;
				power_nl=power_nl&c;
				method=&c;
				lod=&lod;
                iter=&iter;
				keep iter lod method power_nl ; run;

	data power_cc_5; set mean_NL&c;
				power_nl=power_cc;
				method=5;
				lod=&lod;
                iter=&iter;
				keep iter lod method power_nl ; run;


run;


	proc append base=merged_nlmixed_100 data=nlmixed_bias&c  force; run;


	%let c = %eval(&c + 1);
			%let dep = %scan(&multdata, &c);
			%end;
		%mend;
		%nl_ttest(multdata= Ystar2 Ystar3   end);**Macro run statement;


data power_nlmixed_&lod;
merge  power_nlmixed1 power_nlmixed2  power_cc_5; by method;run;
	


	%let d = %eval(&d + 1);
	%let file = %scan(&NLdata, &d);
	%end;
%mend;
%NL(NLdata=   POWER_10018  end);
/**/
/*proc datasets nolist lib=work;*/
/*delete merged_nlmixed_107:  ; quit;*/
/*=================================== POWER of the t-test =============================*/
proc print data=power_nlmixed_12;run;

data power_nlmixed_12;
merge  power_nlmixed1 power_nlmixed2  power_cc_5; by method;run;

/*============================== create merging power data ============================*/

data power_nlmixed_100; merge power_nlmixed_8 power_nlmixed_10 power_nlmixed_12 power_nlmixed_14 power_nlmixed_16 power_nlmixed_18; by lod;run; 

data power_nlmixed_110; merge power_nlmixed_8 power_nlmixed_10 power_nlmixed_12 power_nlmixed_14 power_nlmixed_16 power_nlmixed_18; by lod; run;


 data merged_bias_110; merge nlmixed_bias_110_8 merged_nlmixed_110; by lod; run;

 proc print data=merged_nlmixed_110;run;

data missing_nodif;
input lod percent ;
datalines;
8   8
10  17
12  30 
14  46
16  62
18  75
;
run;

/*******************power NLMixed**********************/
proc sql; create table sasdata.adhoc_nlmixed_110 as 
select merged_bias_110.* , missing_nodif.lod , missing_nodif.percent 
from merged_bias_110, missing_nodif
where merged_bias_110.lod = missing_nodif.lod; quit;
/*==================================================proc format ==================================*/


/****************PROC format****************/
proc format;
    value method_nlmixed
        1="LOD-Nlmixed"
        2="(LOD/2)-NLmixed"
        3="(1/sqrt2)*(LOD)-NLmixed"
        4="0-NLmixed";
	
      
run;
data sasdata.adhoc_nlmixed_105; set sasdata.adhoc_nlmixed_105; 
format method_nl method_nlmixed; run;

/*=============================================================================================*/
/*=========================== plot the adhoc imputation technique =============================*/

title 'Boxplot to display bias under different hybrid imputation techniques over different LOD values and percentage of missing data';
proc sgpanel data= sasdata.adhoc_nlmixed_110;  panelby method_nl /
layout= panel;
format method_nl method_nlmixed.;
   vbox  Bias_nl / category=percent group=LOD  ; 
 refline 0 1 -1 /axis=Y lineattrs=(pattern=dash); 

 colaxis label="Bias under simple combined with nlmixed imputation techniques  at gamma=1.10 ";
 rowaxis label="Range of Bias and its median and IRQ" values=(-4 to 2 by 2);
run; title;



/*=================================== POWER of the t-test =============================*/

proc sql; create table power_nl_100 as 
select power_nl_100.* , missing_nodiffuck.lod , missing_nodiffuck.percent 
from power_nl_100, missing_nodiffuck
where power_nl_100.lod = missing_nodiffuck.lod; quit;

/*================================= POWER OF THE T-TEST FORMAT===============================*/


proc format;
    value method_nlmixed
        1="LOD-Nlmixed"
        2="(LOD/2)-NLmixed"
        3="(1/sqrt2)*(LOD)-NLmixed"
        4="0-NLmixed"
5="complete-cases";
	
      
run;
data sasdata.power_nlmixed_110; set sasdata.power_nlmixed_110; 
format method method_nlmixed; run;


/*================================= PLOT OF POWER T-TEST ===============================*/

  title 'power plot '; 
proc sgplot data=sasdata.power_nlmixed_105 ;
format method method_nlmixed.;

series x=percent y=power_nl / lineattrs=(pattern=solid thickness=2) group=method groupdisplay=cluster clusterwidth=0.5 name="S" markers ;

refline 0.9 /axis=Y lineattrs=(pattern=dash); 
refline 0.1 /axis=Y lineattrs=(pattern=dash); 
yaxis label='power of the t-test under combined adhoc and nlmixed imputation techniques at gamma=1.05 ' values=(0 to 1 by 0.02);
xaxis display=(nolabel);

run;title;



/*============================= Type 1 error data and tabulation ===========================*/
proc print data=sasdata.final_power_100;run;
proc print data=sasdata.power_mixed_100;run;
proc print data=sasdata.power_nlmixed_100;run;
proc print data=sasdata.power_ad_nl_100;run;

data power_nl_100; set power_nl_100; where(method= 7) or (method=1);
if method= 1 then new_method=7; else new_method=7; 
power= power_nl;
drop power_nl method;
run;


data adhoc_100; set sasdata.final_power_100; where (method= 2) OR (method=3);
if method=2 then new_method= 1; else new_method=2; drop method; run;
proc sort data=adhoc_100; by new_method; run;
data mixed_100; set sasdata.final_power_100; where (method=2) OR (method=3); 
if method=2 then new_method=3; else new_method=4; drop method; run;
proc sort data=mixed_100; by new_method;run;
data nlmixed_100; set power_nlmixed_100;  where (method= 1) OR (method=2);
if method=1 then new_method=5; else new_method=6; 
power=power_nl;
drop power_nl method;
run;
proc sort data=nlmixed_100; by new_method; run;
/*========================= merging datasets =======================*/

data sasdata.Type1_error; merge adhoc_100 mixed_100 nlmixed_100 power_nl_100; by new_method; run;

data type1_error; set sasdata.type1_error; type_1_error= power;run;
/*========================= proc univariate ========================*/



/*====================== proc tabulate =============================*/

proc format;
   value method 1='LOD/2'
                2='(1/sqrt2)*LOD'
              
        3="(LOD/2)*Mixed"
        4="(1/sqrt2)*(LOD)*Mixed"
        5="(LOD/2)*Nlmixed"
        6="(1/sqrt2)*(LOD)*Nlmixed"
7= "NlmixedModel";
 
run;
data sasdata.type1_error; set  sasdata.type1_error; 
format new_method method.; run;	






ods rtf file='\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\new_merging.doc' style=statistical;
options nodate pageno=1 linesize=80 pagesize=60;	
proc tabulate data= sasdata.type1_error ;

   class new_method lod percent ;

   var power ;
    table percent,
          new_method*power
          / rts=25 nocontinued;
	
   format new_method method.;

   title 'Energy Expenditures for Each Region';
   title2 '(millions of dollars)';
run; ods rtf close; 