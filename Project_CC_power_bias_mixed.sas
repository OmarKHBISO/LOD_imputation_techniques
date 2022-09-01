
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
					set rep; where (iter=100);
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

ods output "Statistics"=lsmeans_ttest5;
				ods output "T-Tests"=pdiff_ttest5;
				proc ttest data=replacement;
					class trt;
					var AUC5;
					by loop iter;
					where  time = 1;
				run;
				 

				data lsmeans_ttest25 noprint;
					set lsmeans_ttest5;
					where Class= 'Diff (1-2)';
					keep loop iter LowerCLMean UpperCLMean mean StdErr;
					format _all_;
				run;
dm "out;clear;log;clear;";

				data pdiff_ttest25 noprint;
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
					method=5;
					if ProbT <= 0.005 then rejectNull = 1;
					else rejectNull = 0;
					if ProbT = . then rejectNull = .; 	format _all_;
					 




			
/*				proc print data=test_ttest&c;run;*/
                dm "out;clear;log;clear;";
				proc sort data=test_ttest5; by iter loop;run;
/******************things are done here ***************************/



	proc means data=test_ttest5 maxdec=2;where (iter=&iter);
	var rejectNull Mean  Stderr ;
	output out=mean_mixed5 mean=power_mixed5  Mean5 Stderr5 ;
	by iter;
	run; 


	data Bias_th;
set test_ttest5; 

keep iter loop lod method mean rejectNull ;	format _all_; run;





%macro repMixed(multdata);
			%let c=1;
			%let dep = %scan(&multdata, &c, ' ');
			%do %until(&dep = end);
			%put In Mixed: File=&File and Data = &Dep;
dm "out;clear;log;clear;";
				ods output "Estimates"=Pdiff_CS&c;
				ods output "Convergence Status"=ConvergenceCS&c;
				ods output "Fit Statistics"=FitCS&c;
				ods output "lsmeans"=LsMeans&c;
		

				proc mixed data= replacement;
class id trt time;
model &dep= trt time trt*time/solution noint covb; 
repeated time/type=CS subject=id;
/*estimate 'diff_trt' trt 1 -1;*/
estimate 'Y at trt*time'   trt 18 -18 trt*time 1  2  2  2  2  2  2  2  2  1  -1 -2 -2 -2 -2 -2 -2 -2 -2 -1  /divisor=2 cl;

by loop iter;
run;

/*******************************************************************************************
				**************************************************************************
				******************************************************************************
				******************************************************************************************/
dm "out;clear;log;clear;";
				*proc print;*run;
				*ods output close;
				*ods listing;
			
dm "out;clear;log;clear;";
				/*dm "out;clear;log;clear;";*/
				data test_CS&c(rename=(Estimate=Estimate&c StdErr=StdErr&c df=df&c Probt=Probt&c Lower=Lower&c
					Upper=Upper&c rejectNull=rejectNull&c pCS=pCS&c TrueRatio=TrueRatio&c Bias_CS=Bias_CS&c
					widthCL_CS=widthCL_CS&c Coverage_CS=Coverage_CS&c));
					set pdiff_CS&c;
					iter=iter;
					if ProbT <= 0.005 then rejectNull = 1;
					else rejectNull = 0;
					if ProbT = . then rejectNull = .;
					pCS=Probt;
					format _all_; run;

proc sort data=test_CS&c; by iter loop;run;

data SI_mixed&c; merge test_cs&c Bias_th; by loop iter; 

	Bias_cs&c=  Estimate&c - Mean;
run;
				
			
				/*proc print data=test_Cs&c;run;*/
dm "out;clear;log;clear;";

				proc means data=SI_mixed&c maxdec=2;where (iter=&iter);
					var rejectNull&c  Estimate&c Stderr&c Bias_CS&c ;
					output out=mean_Mixed&c mean=power_mixed&c  meanEstimate_CS&c meanStderr_CS&c meanBias_Mixed&c ;
					by iter;
				run;
	
/*************************merging mixed**************************/

	data merged_mixed&c;
				set SI_mixed&c; 
				Bias=Bias_CS&c;
				Method= &c;
				estimate= estimate&c;
				stderr= stderr&c;
			    lod=&lod;
				keep iter loop method bias stderr  estimate lod;
				run;

	data power_m&c; set mean_mixed&c;
				power=power_mixed&c;
				method=&c;
				lod=&lod;
                iter=&iter;
				keep iter lod method power ; run;

	data power_m5; set mean_mixed5;
	power=power_mixed5;
	method=5;
	lod=&lod;
	iter=&iter;
	keep iter lod method power;run;




/*********************merging estimates by mixed model *********************/
	proc append base=merged_mixed&lod data=merged_mixed&c  force; run;
/********merging all LOD configuration using proc append**********************/
		proc append base=merged_mixed_bias data=merged_mixed&lod  force; run;
dm "out;clear;log;clear;";
			%let c = %eval(&c + 1);
			%let dep = %scan(&multdata, &c);
			%end;
		%mend;
		%repMixed(multdata=Ystar1 Ystar2 Ystar3 Ystar4  end);**Macro run statement;





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
/*proc append base=missing_nodif data=missing_nodif_&lod force;*/
/*run;*/


/*data  merged_adhoc_&lod;*/
/*merge merged_ttest1  merged_ttest2 merged_ttest3 merged_ttest4 merged_ttest5 ;by method; run;*/
/**/
data merged_power_&lod;
merge power_m1 power_m2 power_m3 power_m4 power_m5; by method;run;

data merged_mixed&lod;
merge merged_mixed1 merged_mixed2 merged_mixed3 merged_mixed4; by 
method;
run;

	%let d = %eval(&d + 1);
	%let file = %scan(&repdata, &d);
	%end;
%mend;


%replacement(repdata= power_115118 end);



/*proc datasets nolist lib=work;*/
/*delete merged_mixed_bias: ; quit;*/

proc print data= merged_power_18; run;

proc print data= power_m1; run;

proc print data=merged_mixed16;run;







/************** merging by lod value **********************/
data merged_mixed_115; merge merged_mixed8 merged_mixed10   merged_mixed12  merged_mixed14  merged_mixed16  merged_mixed18 ; by lod;run; 

data power_mixed_115; merge merged_power_8  merged_power_10  merged_power_12 merged_power_14 merged_power_16 merged_power_18; by lod; run;
/********** view the results *************/


/** deleting replicated data **********/
proc print data=missing_nodif; run;
/*data missing_nodif; merge missing_nodif_8 missing_nodif_10 missing_nodif_12 missing_nodif_14 missing_nodif_16 missing_nodif_18; by percent; run; */
/*proc sort data=missing_nodif out=missing_nodif nodupkey; where (corr ne 0.3);*/
/*    by lod;*/
/*run;*/
proc print data=missing_nodif; run;
proc contents data=merged_adhoc;run;

/*proc datasets nolist lib=work;*/
/*delete power_adhoc:  missing_nodif:; quit;*/

/***********************************plotting data ***************************/
proc sql; create table sasdata.final_mixed_120 as 
select merged_mixed_120.* , missing_nodif.lod , missing_nodif.percent 
from merged_mixed_120, missing_nodif
where merged_mixed_120.lod = missing_nodif.lod; quit;



/****************PROC format****************/
proc format;
    value method_mixed
        1="LOD*Mixed"
        2="(LOD/2)*Mixed"
        3="(1/sqrt2)*(LOD)*Mixed"
        4="0*Mixed"
		5="complete data";
      
run;
data sasdata.final_mixed_110; set sasdata.final_mixed_110; 
format method method_mixed; run;


title 'Boxplot to display bias under different hybrid imputation techniques over different LOD values and percentage of missing data';
proc sgpanel data= sasdata.final_mixed_120;where (method ne 5);  panelby method /
layout= panel;
format method method_mixed.;
   vbox  Bias / category=percent group=LOD  ; 
 refline 0 1 -1 /axis=Y lineattrs=(pattern=dash); 

 colaxis label="Bias under different PROC Mixed imputation techniques  at gamma=1.2 ";
 rowaxis values=(-4 to 2 by 2);
run; title;



/****************faire avec mean-BIas wiwth std deviation************/

 
 data final_power_120; set sasdata.final_power_120; gamma=1.2; drop percent;run;

  data final_power_115; set sasdata.final_power_115; gamma=1.15; drop percent;run;
   data final_power_110; set sasdata.final_power_110; gamma=1.10; drop percent;run;
    data final_power_105; set sasdata.final_power_105; gamma=1.05; drop percent;run;


	/************* try to add the percent and merge all data to plot  ********************************/
proc sql; create table final_power_120 as 
select final_power_120.* , missing_nodif.lod , missing_nodif.percent 
from final_power_120, missing_nodif
where final_power_120.lod = missing_nodif.lod; quit;

data sasdata.final_adhoc_power; merge final_power_105 final_power_110 final_power_115 final_power_120 ; by gamma; run;
  
/*************** plot mean estimate I guess *******goptions reset=all;run; ploting ******/
/****************************************************************************************/



/****************** ANALYSE THE POWER OF THE IMPUTATION TECHNIQUES *******************/
/*************************************************************************************/
/*************************************************************************************/
/***********************************plotting data for power ***************************/


proc sql; create table sasdata.power_mixed_115 as 
select power_mixed_115.* , missing_nodif.lod , missing_nodif.percent 
from power_mixed_115, missing_nodif
where power_mixed_115.lod = missing_nodif.lod; quit;


proc sort data=sasdata.power_mixed_115 ; by method lod percent;run;

/*proc means data= sasdata.final_power_130 noprint;*/
/*by method percent ;*/
/*var power;*/
/*output out= power mean=mean_power; run;*/
/****************PROC format****************/
proc format;
  value method_mixed
    1="LOD-Mixed"
       2="(LOD/2)-Mixed"
    3="(1/sqrt2)*(LOD)-Mixed"
    4="0-Mixed"
		5="complete data";
      
run;
data sasdata.final_mixed_power; set sasdata.final_mixed_power; 
format method method_mixed; run;

proc print data=sasdata.power_mixed_110;run;

  title 'power plot '; 
proc sgplot data=sasdata.power_mixed_110 ;
format method method_mixed.;

series x=percent y=power / lineattrs=(pattern=solid thickness=2) group=method groupdisplay=cluster clusterwidth=0.5 name="S" markers ;

refline 0.9 /axis=Y lineattrs=(pattern=dash); 
refline 0.5 /axis=Y lineattrs=(pattern=dash); 
yaxis label='power of the t-test under proc mixed imputation techniques at gamma=1.2 ' values=(0 to 1 by 0.02);
xaxis display=(nolabel);

run;title;

/******** import datasets to try to merge and plot in sgpannel ****************/


proc print data="//eu.boehringer.com/BIcorp/stv/ah/RDM/Analyses_Statistiques/stages/2022_Kharriche_Omar/Admin/SAS/power_mixed_120.sas7bdat";

run;


/******* trying to merge datasets to get an unique simple imputation dataset **************/
 
 data final_power_120; set sasdata.final_power_120; gamma=1.2; drop percent;run;

  data final_power_115; set sasdata.final_power_115; gamma=1.15; drop percent;run;
   data final_power_110; set sasdata.final_power_110; gamma=1.10; drop percent;run;
    data final_power_105; set sasdata.final_power_105; gamma=1.05; drop percent;run;


	/************* try to add the percent and merge all data to plot  ********************************/
proc sql; create table final_power_120 as 
select final_power_120.* , missing_nodif.lod , missing_nodif.percent 
from final_power_120, missing_nodif
where final_power_120.lod = missing_nodif.lod; quit;

data sasdata.final_adhoc_power; merge final_power_105 final_power_110 final_power_115 final_power_120 ; by gamma; run;
/********************* trying to merge datasets to get a unique data with different gamma configuration **********/
data power_mixed_120; set power_mixed_120; gamma=1.2;run;

data power_mixed_115; set power_mixed_115; gamma=1.15;run;
data power_mixed_110; set sasdata.power_mixed_110; gamma=1.1;run;
data power_mixed_105; set sasdata.power_mixed_105; gamma=1.05;run;

/*================================================================================================*/
/************* try to add the percent and merge all data to plot  ********************************/
proc sql; create table power_mixed_120 as 
select power_mixed_120.* , missing_nodif.lod , missing_nodif.percent 
from power_mixed_120, missing_nodif
where power_mixed_120.lod = missing_nodif.lod; quit;

/************** plot the power at each size gamma ************************/
data sasdata.final_mixed_power; merge power_mixed_105 power_mixed_110 power_mixed_115  power_mixed_120; by gamma; run;



/************* i have data of new power configuration at 1.1 1.05 1.15 and      and i need the oder configurations to get an sgpannel*/

/*** problem with permission ****/
ods listing style=journal;
title "Power of the t-test";
proc sgpanel data= sasdata.final_adhoc_power; 
styleattrs datasymbols=(circlefilled trianglefilled diamondfilled triangledownfilled circle );
panelby gamma / layout=panel novarname uniscale=column;
series x=percent y=power / group=method lineattrs= (pattern=solid)
markers markerattrs= (color=cx2f2f2f size=8) name='a';
colaxis values=(1 to 5 by 1) integer label='fuck u';
rowaxis offsetmax= .1 label= "values" grid; keylegend 'a' / title="patient:";run; 




title 'Boxplot to display different imputation techniques over different LOD values and percentage of missing data';
proc sgpanel data= sasdata.final_mixed_power;where(method ne .);  

styleattrs datasymbols=(circlefilled trianglefilled diamondfilled triangledownfilled circle );
panelby gamma /layout= panel novarname uniscale=column;
format method method_mixed.;

series x=percent y=power /  group=method  lineattrs = (thickness = 2) markers markerattrs=(symbol=circle size=5 ) name='a';
refline 0 /axis=Y lineattrs=(pattern=dash); 
 colaxis label=  'Sgpanel of power with different mixed imputation techniques vs complete data at increasing size effect gamma and proportion of missing values';
 rowaxis values=(0 to 1 by 0.02) grid;
run; title;



/*********creating an XLS library****************************/
libname sasdata xlsx '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\sasdata.xlsx';

/*=============================================================*/
/*=====================export datasets=========================*/
/*=============================================================*/
proc print data=sasdata.final_adhoc_power (obs=10);run;
	proc export data=sasdata.final_mixed_120 outfile= '//eu.boehringer.com/BIcorp/stv/ah/RDM/Analyses_Statistiques/stages/2022_Kharriche_Omar/Admin/SAS/new_data/sasdata.final_mixed_power.xlsx'
dbms=xlsx; run;

