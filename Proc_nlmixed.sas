
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
		set safe.&file;where (iter=200);


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

	
			%put _USER_;
				%put In NLMixed: File=&File;
		

		
dm "out;clear;";
ods output additionalestimates= Pdiff_NL;
ods output ConvergenceStatus=Convergence_NL;
proc nlmixed data=NL; 

bounds sigsq1 sigsqe >0;
pi= 2*arsin(1);
mu= beta1*trtA + beta2*trtB + beta3*t1 + beta4*t2 + beta5*t3 + beta6*t4 + beta7*t5 + beta8*t6 + beta9*t7 + beta10*t8 + beta11*t9 + beta12*t10 +
beta13*int1 + beta14*int2 + beta15*int3 + beta16*int4 + beta17*int5 + beta18*int6 +beta19*int7 + beta20*int8 + beta21*int9 + beta22*int10 + beta23*int11 + beta24*int12
				+beta25*int13+beta26*int14+beta27*int15+beta28*int16+ beta29*int17+ beta30*int18+ beta31*int19+ beta32*int20+ai;
if conc  ne . then ll = (1/sqrt(2*pi*sigsqe))*exp(-(conc-mu)**2/(2*sigsqe));
if conc  = . then ll = probnorm((conc-mu)/sqrt(sigsqe));
	L = log(ll);
model conc ~ general (L);

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

	data test_NL;
			set pdiff_NL;
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

data nlmixed_bias; merge test_nl test_ttest5; by loop iter ;run;


data nlmixed_bias; set nlmixed_bias;
	Bias_NL = Estimate - mean;
	method_nl =7;

	
run;









		proc print data=nlmixed_bias;run;

		proc means data=nlmixed_bias	 maxdec=2;
		var rejectNull_1 rejectnull Estimate mean Stderr Bias_NL lod method_nl;
	output out=mean_NL mean= power_cc power_nl  meanEstimate_NL meanEstimate_cc meanStderr_NL meanBias_NL lod method_nl;
		by iter;
run;

	data power_nlmixed; set mean_NL;
				power_nl=power_nl;
				method=7;
				lod=&lod;
                iter=&iter;
				keep iter lod method power_nl ; run;

	data power_cc_5; set mean_NL;
				power_nl=power_cc;
				method=5;
				lod=&lod;
                iter=&iter;
				keep iter lod method power_nl ; run;


run;







data power_nl_&lod;
merge  power_nlmixed  power_cc_5; by method;run;
	
	proc append base=merged_nlmixed_100 data=nlmixed_bias  force; run;

	%let d = %eval(&d + 1);
	%let file = %scan(&NLdata, &d);
	%end;
%mend;
%NL(NLdata=   POWER_10018  end);

/*=============================================== POWER NLMIXED ===========================*/

proc print data=power_nlmixed_8;run;

data power_nl_100; merge power_nlmixed_8 power_nl_10 power_nl_12 power_nl_14 power_nl_16 power_nl_18; by lod; run;
