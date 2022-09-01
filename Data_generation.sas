/*options nonotes nodate nosource nosource2 errors=0;*/

proc datasets nolist lib=safe;
delete p: c:  ; quit;
proc datasets nolist lib=work;
delete p_105_8: m_105_8: ; quit;
/*================================================================================================================================*/
/*============================================= SIMULATION OF LONGITUDINAL DATA ==================================================*/
/*================================================================================================================================*/

libname sasdata '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS';
libname safe '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\LOD_Methods_code\new_datasets';
%GLOBAL n iter multfactor factor ddfm letmesee rho LOD;


%let n = 40;
%let iter =1;
%let multfactor= 100;
%let rho= 1;
%PUT _USER_;
*options MPRINT;
%macro gibaldi_cv();
%do LOD=8 %to 18 %by 2;   /*%do LOD=16 %to 22 %by 2;*/
      /*%do rho=3 %to 7 %by 4;*/
/*dm "out;clear;log;clear;";*/
	%let c=1;
	%let dep = %scan(&multfactor, &c, ' ');
	/*%do %until(&dep = 500);*/
/*	dm "out;clear;log;clear;";*/
		%let factor = %scan(&multfactor, &c, ' ');
		%let letmesee=%sysevalf(&factor/100);
		%do loop = 1 %to &iter;

			proc iml;
				time=j(1,10,.);
				time = {1 2 3 4 5 6 7 8  9  10};
				p = ncol(time);
				LOD=&LOD;
				n=&n;

				call symput('p',left(char(p)));call symput('n',left(char(n)));
				muH0=j(p,1,.);muHA=j(p,1,.);
				z=j(p,1,.);miss=j(n,p,0);
				R=j(p,p,.);H=j(p,p,.);
						H={ 0    1    1    1    1    1    1    1   1   1,
				    1    0    1    1    1    1    1    1   1   1,
				    1    1    0    1    1    1    1    1   1   1,
				    1    1    1    0    1    1    1    1   1   1,
				    1    1    1    1    0    1    1    1   1   1,
				    1    1    1    1    1    0    1    1   1   1,
                    1    1    1    1    1    1    0    1   1   1,
                    1    1    1    1    1    1    1    0   1   1,
                    1    1    1    1    1    1    1    1   0   1,
                    1    1    1    1    1    1    1    1   1   0};
				%let rhodec=%sysevalf(&rho/10);
				R=&rhodec**H;

				%let BprobDec=%sysevalf(&LOD/10);
				BprobDec=&BprobDec;
	/*muH0={1.2,1.9,2.1,2.3,2.5,3,2.5,1.5,1.5,1.5};*/
		muH0={1,1.2,1.4,1.6,2,1.8,1.6,1.4,1.2,1.1};

			muH1={1,1.5,2,2.2,2.5,3,2.5,2,1.7,1.5};
			/*muH1={1.1,1.9,2.1,2.4,2.6,3,2.5,1.5,1.6,1.5};*/
				muHA=(muH0*&letmesee);
				sigma=j(p,1,0.361);** so variance = 0.0961! 0.1; **To be used when variance = constant;

				diag_sig = diag( sigma );
				DRD = diag_sig * R * diag_sig`;
				U = half(DRD);

				loop=&loop;iter=&iter;
				%do k = 1 %to (&n/2); 
					z = rannor(J(p,1,0));
				/*	y = muH0 + ((U` * z)-(U` * z));*/
					y = muH0 + (U` * z);
					yprime = y`;
					yall = yall // yprime;
				%end;
				%do k = ((&n/2)+1) %to &n; 
					z = rannor(J(p,1,0));
		/*	y = exp(muHA + ((U` * z)-(U` * z)));*/
				y = muHA + (U` * z);
					yprime = y`;
					yall = yall // yprime;
				%end;

		varnames = {y1 y2 y3 y4 y5 y6 y7 y8  y9 y10};
				create my_MVN from yall (|colname = varnames|);
				append from yall;
				create time from time;
				append from time;
			quit;
			*proc print data=my_MVN;
			*run;
			*proc print data=time;
			*run;

			data mvn ;
				if _n_=1 then set time;
				merge my_MVN;
				id=_n_;
				loop=&loop; *********************************************;
				iter=&iter; *********************************************;
				gamma=&letmesee; ****************************************;
				corr=&rhodec; *******************************************;
				LOD=&LOD; *******************************************;
				t=0;Y=Y1;time=Col1;output;
				t=1;Y=Y2;time=Col2;output;
				t=2;Y=Y3;time=Col3;output;
				t=3;Y=Y4;;time=Col4;output;
				t=4;Y=Y5;time=Col5;output;
				t=5;Y=Y6;time=Col6;output;
				t=6;Y=Y7;time=Col7;output;
				t=7;Y=Y8;time=Col8;output;
				t=8;Y=Y9;time=Col9;output;
				t=9;Y=Y10;time=Col10;output;
				keep id t time y LOD loop gamma corr iter ;
			run;

			data crossover ;
				set mvn;
				if (id <= 20) then do;
					conc=abs(Y); if abs(Y) < (LOD/10) then conc = .;
					trt='A';seq=1;			
				end;
				if (id > 20) then do;
					conc=abs(Y);if abs(Y) < (LOD/10) then conc = .;
					trt='B';seq=2;
				end;
			run;
			/* creation de la data complete cases */
				data crossover2 ;
				set mvn;
				if (id <=20) then do;
					conc=abs(Y);
					trt='A';seq=1;			
				end;
				if (id > 20) then do;
					conc=abs(Y);
					trt='B';seq=2;
				end;
			run;
			proc sort data=crossover;
				by loop trt;
			run;
/* sort CC data */
	proc sort data=crossover2;
				by loop trt;
			run;
			data AUC;
				set crossover;
	           if (t=0|t=9) then do;
					if conc=. then conc=0;
				end;
				if (conc=.) then delete;
				drop LagTime LagValue;
				LagTime = LAG(time);
				LagValue = LAG(conc);

				if t = 0 then do;
					LagTime = 0;
					LagValue = 0;
					SumTrapezoid=0; **My addition: Does AUC for 'by id';
					Trapezoid = .;*My addition: Make sure that trapezoid = 0 at first time point;
					SumTrapezoid + Trapezoid;
				end;
				else do;
					Trapezoid = (time-LagTime)*(conc+LagValue)/2; 
					SumTrapezoid + Trapezoid;
				end;
			run;

		/* AUC calcul for CC data */
				data AUC2;
				set crossover2;
	           if (t=0|t=9) then do;
					if conc=. then conc=0;
				end;
			
				LagTime = LAG(time);
				LagValue = LAG(conc);

				if t = 0 then do;
					LagTime = 0;
					LagValue = 0;
					SumTrapezoid2=0; **My addition: Does AUC for 'by id';
					Trapezoid2 = .;*My addition: Make sure that trapezoid = 0 at first time point;
					SumTrapezoid2 + Trapezoid2;
				end;
				else do;
					Trapezoid2 = (time-LagTime)*(conc+LagValue)/2; 
					SumTrapezoid2 + Trapezoid2;
				end;
			run;
			/*proc print data=crossover;run;
			proc print data=AUC;run;*/
            data AUC_global; 
			merge AUC AUC2;  by id loop trt iter LOD corr gamma; run;

			data AUCBySubject;
				set AUC_global;
				if time = 10;
				AUC=SumTrapezoid2;
				keep id AUC trt loop gamma corr iter LOD;
			run;

			proc sort data=AUCBySubject; by id;run;
			proc sort data=crossover;by id;run;
			
			/*proc print data=crossover;run;*/
			/*proc print data=AUCBySubject;run;*/
			data all  ;
				merge crossover AUCBySubject;
				by id;
			run;
			/*proc print data=all;run;*/
			proc append base=safe.power_&Factor&rho&LOD data=all force;
			run;
		
			/*proc print data=safe.power_100319;run;*/
/*			dm "out;clear;log;clear;";*/
		%end;

		/*Analysis to be entered here*/

	    %let c = %eval(&c + 1);
	    %let dep = %scan(&multfactor, &c);
	%end;

%mend;
%gibaldi_cv(); /*100 102 105 107 110 500*/ **Macro run statement;

%macro obtainDataSets();
	ods output members=DataSetListTemp(Keep=name LastModified);
	proc datasets mt=data library=safe;
	run;
	quit;

	proc print data=DataSetList noobs;
	run;



%mend;
%obtainDataSets();

  /*************/
/**********************************************/
/**********************************************/
/**********************************************/
/* plot longitudinal mesures Y by id and time */
/**********************************************/
/**********************************************/
/**********************************************/
	proc export data=safe.power_115118 outfile= '//eu.boehringer.com/BIcorp/stv/ah/RDM/Analyses_Statistiques/stages/2022_Kharriche_Omar/Admin/SAS/new_data/safe.power_115118.xlsx'
dbms=xlsx; run;












/*proc print data= work.power_58315; run;*/

title "longitudinal distribution of viral shedding by treatment group";
proc sgplot data=all; loess x=t y=Y/group= trt clm clmtransparency=0.6; run;title;

data all; set all; id_trt= id||trt; run;

title "longitudinal distribution of viral shedding by treatment group";
proc sgpanel data=all;
	styleattrs datalinepatterns=(solid solid);
	panelby trt/ rows=2 novarname;
	format id_trt $Group.;
	series x=Time y=conc/ group=id_trt transparency=0 markers lineattrs = (thickness = 1 pattern = dot);
/*	series x=Time y=Y/  group=id_trt lineattrs = (thickness = 1 pattern = solid);*/

	rowaxis label="Response Y (viral shedding)";
	
run;title;
/*************************plot data under each LOD value in group A **************************/
title "longitudinal distribution of viral shedding under LOD=1.6 for Treatment A";
proc sgplot data=all; where trt="A";
series x= time y=Conc/ group=id_trt transparency=0 markers lineattrs = (thickness = 1);
yaxis min=0 max=4;
refline 1.6 / axis=y lineattrs=(thickness=2 color=darkblue pattern=dash );
run;title;
/*************************plot data under each LOD value in group B **************************/
title "longitudinal distribution of viral shedding under LOD=1.6 for Treatment B";
proc sgplot data=all; where trt="B";
series x= time y=Conc/ group=id_trt transparency=0 markers lineattrs = (thickness = 1);
yaxis min=0 max=4;
refline 1.6 / axis=y lineattrs=(thickness=2 color=darkblue pattern=dash );
run;title;


proc export data=safe.power_13038 outfile = '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\power_13038.xls' dbms=xls replace;run;
proc export data=safe.power_130310 outfile = '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\power_130310.xls' dbms=xls replace;run;
proc export data=safe.power_130312 outfile = '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\power_130312.xls' dbms=xls replace;run;
proc export data=safe.power_130314 outfile = '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\power_130314.xls' dbms=xls replace;run;
proc export data=safe.power_130316 outfile = '\\eu.boehringer.com\BIcorp\stv\ah\RDM\Analyses_Statistiques\stages\2022_Kharriche_Omar\Admin\SAS\power_130316.xls' dbms=xls replace;run;

