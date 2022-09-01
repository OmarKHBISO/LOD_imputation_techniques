/*proc datasets nolist lib=work;*/
/*delete missing_nodif: ; quit;*/
/**/
/*proc datasets lib=sasdata KILL;*/
/*quit;*/

/*options nonotes nodate nosource nosource2 errors=0;*/
/*dm "out;clear;log;clear;";*/
/*dm "out;clear;log;clear;";*/
libname sasdata '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS';
libname safe '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\LOD_Methods_code\new_datasets';
%GLOBAL n iter multfactor factor ddfm letmesee rho LOD;

%PUT _USER_;


%macro replacement(repdata);
	%let d=1;
	%let file = %scan(&repdata, &d, ' ');
	%do %until(&file = end);




	/*Analysis begins here*/
	
dm "out;clear;log;clear;";
		data rep noprint;
			set safe.&file;
			/*set safe.POWER_125330;*/
		
		/*	if conc = . then Ystar1 = .; else Ystar1 = conc;*/
			if conc = . then Ystar1 = (1*(LOD/10)); else Ystar1 = conc;
			if conc = . then Ystar2 = ((1/2)*(LOD/10)); else Ystar2 = conc;
			if conc = . then Ystar3 = ((1/sqrt(2))*(LOD/10)); else Ystar3 = conc;
			if conc = . then Ystar4 = ((0)*(LOD/10)); else Ystar4 = conc;
		   	Ystar5= Y;
			*if loop = 5;
			*where loop in (1:1000);
		run;
	
dm "out;clear;log;clear;";
		%macro AUC(multdata);
			%let c=1;
			%let dep = %scan(&multdata, &c, ' ');
			%do %until(&dep = end);

 
				data AUCRep&c noprint;
					set rep; where (iter=200);
					if (t=0|t=9) then do;
						if &dep=. then &dep=0;
					end;
					if (&dep=.) then delete;
				/*	drop LagTime LagValue;*/
					LagTime = LAG(time);
					LagValue = LAG(&dep);

					if t = 0 then do;
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
dm "out;clear;log;clear;";
			%let c = %eval(&c + 1);
		    %let dep = %scan(&multdata, &c);
			%end;
			/*proc print data=AUCRep1;run;*/
			/*proc print data=AUCRep2;run;*/
			/*proc print data=AUCRep3;run;*/
			/*proc print data=AUCRep4;run;*/
			/*proc print data=AUCRep5;run;*/
dm "out;clear;log;clear;";
		    proc sort data=AUCRep1;by loop id;run;
			proc sort data=AUCRep2;by loop id;run;
			proc sort data=AUCRep3;by loop id;run;
			proc sort data=AUCRep4;by loop id;run;
		    proc sort data=AUCRep5;by loop id;run;
		
					
			data AUCRep noprint;
				merge  AUCRep1 AUCRep2 AUCRep3 AUCRep4 AUCRep5;
				by loop id;
			run;
dm "out;clear;log;clear;";
			/*proc print data=AUCRep;run;*/
			data AUCBySubjectRep;
				set AUCRep;
				if time = 10;
				/* here i need to calculate the AUC_CC from the value of Y */
	

AUC1=SumTrapezoid1;AUC2=SumTrapezoid2;AUC3=SumTrapezoid3;AUC4=SumTrapezoid4; AUC5=SumTrapezoid5;
				keep id AUC1 AUC2 AUC3 AUC4  AUC5 trt loop gamma corr iter LOD;
			run;
dm "out;clear;log;clear;";
		
			/*proc print data=AUCRep;run;*/
			proc sort data=AUCBySubjectRep; by loop id;run;
			proc sort data=AUCRep;by loop id;run;
dm "out;clear;log;clear;";
			data replacement;
				merge AUCRep AUCBySubjectRep;
				by loop id;
			run;
			/*proc print data=replacement;run;*/
dm "out;clear;log;clear;";
		%mend;
		%AUC(multdata=   Ystar1 Ystar2 Ystar3 Ystar4  Ystar5 end); **Macro run statement;
		**Check if merge above is right;
dm "out;clear;log;clear;";
		data _null_;
			set replacement;
			call symput('gamma',left(gamma));
			call symput('LOD',left(LOD));
			call symput('corr',left(corr));
			call symput('iter',left(iter));
		run;
		%put _USER_;
dm "out;clear;log;clear;";
		proc sort data=replacement;
			by loop iter;
		run;
		/*proc print data=replacement;run;*/
		/*proc means data=replacement;
			class time;
			var conc;
			where loop = 5;
		run;*/
dm "out;clear;log;clear;";
		%macro repTTest(multdata);
			%let c=1;
			%let dep = %scan(&multdata, &c, ' ');
			%do %until(&dep = end);
dm "out;clear;log;clear;";
				ods output "Statistics"=lsmeans_ttest&c;
				ods output "T-Tests"=pdiff_ttest&c;
				proc ttest data=replacement;
					class trt;
					var &dep;
					by loop iter;
					where  time = 1;
				run;
				 
dm "out;clear;log;clear;";
				data lsmeans_ttest2&c noprint;
					set lsmeans_ttest&c;
					where Class= 'Diff (1-2)';
					keep loop iter LowerCLMean UpperCLMean mean StdErr;
					format _all_;
				run;
dm "out;clear;log;clear;";

				data pdiff_ttest2&c noprint;
					set pdiff_ttest&c;
					where Method='Satterthwaite';
					keep loop iter df Probt;
					format _all_;
				run;
dm "out;clear;log;clear;";
				/*dm "out;clear;log;clear;";*/
				*ods listing close;
				data test_ttest&c(rename=(df=df&c Probt=Probt&c Mean=Mean&c LowerCLMean=LowerCLMean&c UpperCLMean = UpperCLMean&c
					StdErr=StdErr&c rejectNull=rejectNull&c TrueRatio=TrueRatio&c Bias=Bias&c
					widthCL_ttest=widthCL_ttest&c Coverage_ttest=Coverage_ttest&c));
					merge pdiff_ttest2&c lsmeans_ttest2&c;
					by loop ;
					iter=&iter;
					if ProbT <= 0.005 then rejectNull = 1;
					else rejectNull = 0;
					if ProbT = . then rejectNull = .; 	format _all_;
					 




			
/*				proc print data=test_ttest&c;run;*/
                dm "out;clear;log;clear;";
				proc sort data=test_ttest&c; by iter loop;run;
/******************things are done here ***************************/
data Bias_th;
set test_ttest5; 
keep iter loop mean5 rejectNull5;	format _all_; run;

proc sort data=Bias_th; by iter loop;run;
				data SI_ttest&c; merge test_ttest&c Bias_th; by loop iter; 
	Bias&c=  Mean&c - Mean5; 
run;
			/*	data SI_ttest1; merge test_ttest1 test_ttest5; by loop iter; 
			Bias1=  Mean1 - Mean5; run;

				data SI_ttest2; merge test_ttest2 test_ttest5; by loop iter; 
			Bias2=  Mean2 - Mean5; run;

				data SI_ttest3; merge test_ttest3 test_ttest5; by loop iter; 
			Bias3=  Mean3 - Mean5; run;
				data SI_ttest4; merge test_ttest4 test_ttest5; by loop iter; 
			Bias4=  Mean4 - Mean5; run;

                 data SI_ttest5; set test_ttest5; Bias5=0;run;*/

				proc means data=SI_ttest&c maxdec=2;where (iter=&iter);

					var rejectNull&c Mean&c Stderr&c Bias&c ;
					output out=mean_ttest&c mean=type1error&c meanEstimate_ttest&c meanStderr_ttest&c meanBias&c ;
					by iter;
				run;
/*********************merging_ttest**************************/
						data merged_ttest&c;
				set SI_ttest&c; 
				Bias=Bias&c;
				Stderr= Stderr&c;
				Method= &c;
			
					stderr= stderr&c;
				lowerCL= lowerCLMEAN&c;
				upperCL= upperCLMEAN&c;
				Mean= Mean&c;
			    lod=&lod;
				keep loop iter method bias  stderr lowercl uppercl mean lod;
				run;
	/*********************power data********************************/
	data power_adhoc&c; set mean_ttest&c;
				power=Type1error&c;
				method=&c;
				lod=&lod;
                iter=&iter;
				keep iter lod method power ; run;
dm "out;clear;log;clear;";
/********try to merge here **********************/

/*	proc append base=sasdata.merged_adhoc_nodif&lod data=merged_ttest&c force;run;*/
/********merging all LOD configuration using proc append**********************/
/*		proc append base=merged_adhoc_nodif data=sasdata.merged_adhoc_nodif&lod  force; run;*/

/*****************merging power ************************/
	
    proc append base=power_adhoc data=merged_power_&lod force;run;
			%let c = %eval(&c + 1);
			%let dep = %scan(&multdata, &c);
			%end;
		%mend;
		%repTTest(multdata=AUC1 AUC2 AUC3 AUC4  AUC5 end); **Macro run statement;
/*******************proportion of censoring *******************/
proc freq data=safe.&file;where (iter=100); table conc/missing out= missing_out	maxlevels=1 ;
run;
/**/
data missing_nodif_&LOD ;
set missing_out (obs=1);
corr=&corr;
LOD=&LOD;
run;
proc print data=missing_nodif_&lod;run;
/****************  Appending SAS missing proportion data **************/
	proc append base=missing_nodif data=missing_nodif_&lod force;
run;


data  merged_adhoc_&lod;
merge merged_ttest1  merged_ttest2 merged_ttest3 merged_ttest4 merged_ttest5 ;by method; run;

data merged_power_&lod;
merge power_adhoc1 power_adhoc2 power_adhoc3 power_adhoc4 power_adhoc5; by method;run;

	%let d = %eval(&d + 1);
	%let file = %scan(&repdata, &d);
	%end;
%mend;


%replacement(repdata= power_120118 end);


/*===================================================================================================*/
/*===================================================================================================*/
/*===================================================================================================*/
/*===================================================================================================*/

proc print data=missing_nodif; run;
proc print data=merged_power_16;run;
proc print data=merged_adhoc_8;run;

/*====================================merging by lod value=============================================*/
data merged_adhoc_120; merge merged_adhoc_8 merged_adhoc_10  merged_adhoc_12 merged_adhoc_14 merged_adhoc_16 merged_adhoc_18; by lod;run; 

data merged_power_120; merge merged_power_8  merged_power_10  merged_power_12 merged_power_14 merged_power_16 merged_power_18; by lod; run;
/********** view the results *************/


/** deleting replicated data **********/
proc sort data=missing_nodif out=missing_nodif nodupkey; where (corr ne 0.3);
   by lod;
run;
/*proc print data=missing_nodif; run;*/
/*proc contents data=merged_adhoc;run;*/
/**/
/*proc datasets nolist lib=work;*/
/*delete power_adhoc:  missing_nodif:; quit;*/

/***********************************plotting data ***************************/
proc sql; create table sasdata.final_adhoc_120 as 
select merged_adhoc_120.* , missing_nodif.lod , missing_nodif.percent 
from merged_adhoc_120, missing_nodif
where merged_adhoc_120.lod = missing_nodif.lod; quit;



/****************PROC format****************/
proc format;
    value method
        1="LOD"
        2="LOD/2"
        3="(1/sqrt2)*(LOD)"
        4="0"
		5="complete data";
      
run;
data sasdata.final_adhoc_115; set sasdata.final_adhoc_115; 
format method method; run;

title 'Boxplot to display bias under different simple imputation techniques over different LOD values and percentage of missing data';
proc sgpanel data= sasdata.final_adhoc_120;where (method ne 5);  panelby method /
layout= panel;
format method method.;
   vbox  Bias / category=percent group=LOD  ; 
 refline 0 1 -1 /axis=Y lineattrs=(pattern=dash); 

 colaxis label="Bias under different imputation techniques at gamma=1.20";
run; title;



/****************faire avec mean-BIas wiwth std deviation************/
proc print data=sasdata.final_adhoc_105(obs=10);run;

  
/*************** plot mean estimate I guess *******goptions reset=all;run; ploting ******/
/****************************************************************************************/



/****************** ANALYSE THE POWER OF THE IMPUTATION TECHNIQUES *******************/
/*************************************************************************************/
/*************************************************************************************/
/***********************************plotting data for power ***************************/


proc sql; create table sasdata.final_power_120 as 
select merged_power_120.* , missing_nodif.lod , missing_nodif.percent 
from merged_power_120, missing_nodif
where merged_power_120.lod = missing_nodif.lod; quit;


proc sort data=sasdata.final_power_120; by method lod percent;run;

/*proc means data= sasdata.final_power_130 noprint;*/
/*by method percent ;*/
/*var power;*/
/*output out= power mean=mean_power; run;*/
/*/****************PROC format****************/*/
/*proc format;*/
/*    value method*/
/*        1="LOD"*/
/*        2="LOD/2"*/
/*        3="(1/sqrt2)*(LOD)"*/
/*        4="0"*/
/*		5="complete data";*/
/*      */
/*run;*/


/*=================== plot the power at each size gamma ==========================*/

  title 'power plot '; 
proc sgplot data=sasdata.final_power_120;
format method method.;

series x=percent y=power / lineattrs=(pattern=solid thickness=2) group=method groupdisplay=cluster clusterwidth=0.5 name="S" markers ;

refline 0.9 /axis=Y lineattrs=(pattern=dash); 
refline 1 /axis=Y lineattrs=(pattern=dash); 
yaxis label='power of the t-test under effect size gamma=1.20 by different simple imputation techniques' values=(0 to 1 by 0.02);
xaxis display=(nolabel);

run;title;

/************* i have data of new power configuration at 1.25 and 1.3 1.2 1.15 1.10 and i need the oder configurations to get an sgpannel*/




/*********creating an XLS library****************************/
libname sasdata xlsx '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\sasdata.xlsx';
proc export data=sasdata.final_power_120 
outfile='\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\new_data\final_adhoc_120.xlsx' replace ;
run;

