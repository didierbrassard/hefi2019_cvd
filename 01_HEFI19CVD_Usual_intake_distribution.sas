 /*************************************************************************/
 /*                                                                       */
 /*               HEFI-2019 and CVD risk in the UK Biobank                */
 /*                                                                       */
 /*                  Brassard et al. Am J Clin Nutr 2022                  */
 /*             Code 1: Usual intake distribution estimation              */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 1                               */
 /*                               05AUG2022                               */
 /*                                                                       */
 /*************************************************************************/

/* Indicate location of files (case sensitive) */
	%global path suffix;
	%let path=C:/Users/bradid01/Documents/HEFI19_OUTCOME/;

 /*************************************************************************/
 /*                                                                       */
 /* GENERAL: SET LIBRARIES AND INCLUDE MACROS                             */
 /*                                                                       */
 /*************************************************************************/

/* Define suffix for current model - used for file indexing within the NCI method */
	%let suffix=_distrib ;

/* National Cancer Institute (NCI) macros */
	/* boxcox svy */ %include "&path./Macros/boxcox_survey.macro.v1.2.sas";
	/* std boxcox */ %include "&path./Macros/std_cov_boxcox24hr_conday_minamt_macro_v2.0.sas";
	/* multi MCMC */ %include "&path./Macros/multivar_mcmc_macro_v2.1.sas";
	/* multi dist */ %include "&path./Macros/multivar_distrib_macro_v2.1.sas";
	/* percentiles*/ %include "&path./Macros/percentiles_survey.macro.v1.1.sas";
	
	/* Available at: https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0 */

/* HEFI-2019 scoring algorithm macro */
	/* hefi-2019 */ %include "&path./Macros/hefi2019.scoring.macro.sas";
	
	/* Available at: https://github.com/didierbrassard/hefi2019/ */
	

/* create or indicate library folder */
	options dlcreatedir;
	libname temp "&path./Temp/";
	libname fmtdata "&path./Fmtdata/";
	
/* auto assign proper libraries */
	%macro mcmclib(suffix);
	options dlcreatedir;	
	libname nci "&path./NCI/";                             /* create NCI folder, if need be */
	libname mcmc "&path./NCI/MCMC&suffix./";               /* create MCMC folder, if need be */
	libname chain "&path./NCI/MCMC&suffix./Markovchains/"; /* folder for trace plots */
	libname reslib "&path./NCI/MCMC&suffix./Results/";     /* results */
	libname baselib "&path./NCI/MCMC&suffix./Model/";      /* base estimates */
	libname bootlib "&path./NCI/MCMC&suffix./Bootlib/";    /* bootstrap replicates */
	%mend;
	
	%mcmclib(suffix=&suffix.);

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: PREPARE DATA SETS NEEDED FOR MULTIVARIATE METHOD        */
 /*                                                                       */
 /*************************************************************************/

/* Prepare an input data of eligible sample, dietary intakes and covariates */
	data preNCI(drop=valid24hr included rename=(nomat= subjectid));
		merge fmtdata.screening(keep=nomat included ageat24hr event_primary event_date_primary cens_death rename=(ageat24hr=age))
			  fmtdata.diet24hr(in=foods drop=mean:)
			  fmtdata.date24hr(keep=nomat r24_date: )
			  fmtdata.valid24hr(in=valid keep=nomat valid24hr where=(valid24hr=1))
			  fmtdata.urine_na_24h(keep=nomat urine_na_24h)
			  fmtdata.zdata(keep=nomat sex missing_top2 rename=(sex=_sex))
			  ;
		by nomat ;
	
	* Make numerical sex variable;
		if _sex = "Female" then sex=1;
		if _sex = "Male" then sex=2;
	
	* Categorize age ;
		if age < 51 then agec=1;
		else if 51 <= age < 71 then agec=2;
			else if age >=71 then agec=3;
			
		* dummy variables ;
		agec_2 = 0;
		agec_3 = 0;
		if agec=2 then agec_2 = 1;
		if agec=3 then agec_3 = 1;

	* Output day of the week (for weekend adjustment/balancing);
		
		/* note: SAS weekday are coded as 1=Sunday, 2=Monday, . . . , 7=Saturday */
	
		array r24_wkd(5);
		array _date_(*) r24_date1-r24_date5;
		do i=1 to 5;
			if missing(_date_(i)) then r24_wkd(i)=.;
			else if weekday(_date_(i)) in (1 6 7) then r24_wkd(i)=1;
			else r24_wkd(i)=0;
		end;
	
	* 24hr must be valid and participants part of the eligible sample ;
	if valid & included=1 & missing_top2=0 then output;
	run;

	/* note: n rows = 136,698 */

	proc freq data=_preNCI;
	title1 "NCI Descriptive. Covariate for descriptive model";
	table sex agec event_primary cens_death agec*sex;
	run;
	title1;

/* Output data with 1-row per subject for current replicate, based on preNCI data */
	data subj1rec;
	set preNCI ;
	* constant1 - used to have an intercept in NCI multivariate models;
	constant1=1;
	* replicate - explicit indication that current data is for original sample;
	replicate=0;
	run;

	proc sort data=subj1rec;
	by subjectid;
	run;

/* Make data with m-rows per subject where m = number of 24-h recall (needed FOR NCI)*/
	data datamrec (drop=vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 sugars_free1 
			energy1 plantbev1 pfpb1 rg1 wg1 ssbs1 vf2 pfab2 other2 water2 milk2 mufa2 
			pufa2 sfa2 sugars_free2 energy2 plantbev2 pfpb2 rg2 wg2 ssbs2 vf3 pfab3 
			other3 water3 milk3 mufa3 pufa3 sfa3 sugars_free3 energy3 plantbev3 pfpb3 rg3 
			wg3 ssbs3 vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 sugars_free4 energy4 
			plantbev4 pfpb4 rg4 wg4 ssbs4 vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 
			sugars_free5 energy5 plantbev5 pfpb5 rg5 wg5 ssbs5 r24_wkd1-r24_wkd5 ith);
		* drop recall-specific names ;
		retain replicate subjectid r24_no;
		set subj1rec;
		
		* generic variable names ;
		array foodlist0 (*) vf pfab other water milk mufa pufa sfa sugars_free energy 
			plantbev pfpb rg wg ssbs;
			
		* assign generic variable names to recall specific variable names ;
		
		* Recall #1: output data if there are no missing values  ;
		if nmiss(of vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 sugars_free1 
			energy1 pfpb1 rg1 wg1 ssbs1)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr1 (*) vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 
					sugars_free1 energy1 plantbev1 pfpb1 rg1 wg1 ssbs1;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr1(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=1;
				* weekend indicator ;
				r24_wkd = r24_wkd1;
				* writes observation for current recall ;
				output;
			end;
	
		* Recall #2: output data if there are no missing values  ;
		if nmiss(of vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 sugars_free2 
			energy2 pfpb2 rg2 wg2 ssbs2)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr2 (*) vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 
					sugars_free2 energy2 plantbev2 pfpb2 rg2 wg2 ssbs2;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr2(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=2;
				* weekend indicator ;
				r24_wkd = r24_wkd2;
				* writes observation for current recall ;
				output;
			end;
			
		* Recall #3: output data if there are no missing values  ;
		if nmiss(of vf3 pfab3 other3 water3 milk3 mufa3 pufa3 sfa3 sugars_free3 
			energy3 pfpb3 rg3 wg3 ssbs3)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr3 (*) vf3 pfab3 other3 water3 milk3 mufa3 pufa3 sfa3 
					sugars_free3 energy3 plantbev3 pfpb3 rg3 wg3 ssbs3;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr3(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=3;
				* weekend indicator ;
				r24_wkd = r24_wkd3;
				* writes observation for current recall ;
				output;
			end;
	
		* Recall #4: output data if there are no missing values  ;
		if nmiss(of vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 sugars_free4 
			energy4 pfpb4 rg4 wg4 ssbs4)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr4 (*) vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 
					sugars_free4 energy4 plantbev4 pfpb4 rg4 wg4 ssbs4;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr4(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=4;
				* weekend indicator ;
				r24_wkd = r24_wkd4;
				* writes observation for current recall ;
				output;
			end;
			
		* Recall #5: output data if there are no missing values  ;
		if nmiss(of vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 sugars_free5 
			energy5 pfpb5 rg5 wg5 ssbs5)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr5 (*) vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 
					sugars_free5 energy5 plantbev5 pfpb5 rg5 wg5 ssbs5;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr5(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=5;
				* weekend indicator ;
				r24_wkd = r24_wkd5;
				* writes observation for current recall ;
				output;
			end;
	run;
	
/* Create new dummy variables for sequence of 24-h recalls */
	data datamrec;
	set datamrec;
	* updated seq indicator ;
	if r24_no = 2 then seq2=1;
	else seq2=0;
	
	if r24_no = 3 then seq3=1;
	else seq3=0;
	
	if r24_no = 4 then seq4=1;
	else seq4=0;
	
	if r24_no = 5 then seq5=1;
	else seq5=0;
	run;

	proc sort data=datamrec;
		by subjectid ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*  PRE NCI STEP: KEEP ONLY PARTICIPANTS IN CURRENT STRATUM (SEX=MALES)  */
 /*                                                                       */
 /*************************************************************************/

	data _tempstratum1 ;
		set datamrec ;
	* keep only data for current stratum (i.e., males) ;
		if (sex = 1) then output;
	run;
	
	proc sort ;
		by subjectid ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* CURRENT STRATUM = MALES ONLY                                          */
 /* PRE NCI STEP: GET MIN AMOUNT AND FIND BEST BOX-COX TRANSFORMATIONS    */
 /*                                                                       */
 /*************************************************************************/

	/* Output data from the first 24-h dietary recall completed only */
	proc sort data = _tempstratum1 nodupkey out=inboxcox;
	by subjectid ;
	run;
	
	/* make a list of covariates for consistency */
	%let covars  = seq2 seq3 seq4 seq5 r24_wkd urine_na_24h agec_2 agec_3 event_primary cens_death ;
	
	/* ensure that the macro variable <best_lambda> is available outside the <boxcox_survey> macro */
	%global best_lambda ;

	/********************************************************************/
	/*** Dietary constituents #1 - plantbev: min. amount and lambda   ***/
	/********************************************************************/

	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var plantbev ;
	where plantbev > 0;
	output out=_min1(keep=minamount) min=minamount ;
	run;
	
	data _min1;
	set _min1;
	minamount=minamount/2;
	tran_paramindex=1;
	length varname $ 32;
	varname="plantbev";
	run;

	data _temp1;
	set inboxcox(keep=subjectid plantbev &covars);
	if (plantbev > 0) then output;
	run;

	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp1, 
	 subject = subjectid, 
	 var     = plantbev,  
	 covars  = &covars. ,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;

	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x1;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	* tran_lambda=0.29; * = value of macro variable for replicate=0, stratum=males ;
	tran_paramindex=1;
	varname="plantbev";
	run;

	/* note: steps 1-3 are repeated for all other dietary constituents of the HEFI-2019 below */

	/********************************************************************/
	/*** Dietary constituents #2 - pfpb: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var pfpb ;
	where pfpb > 0;
	output out=_min2(keep=minamount) min=minamount ;
	run;
	
	data _min2;
	set _min2;
	minamount=minamount/2;
	tran_paramindex=2;
	length varname $ 32;
	varname="pfpb";
	run;
	
	data _temp2;
	set inboxcox(keep=subjectid pfpb &covars. );
	if (pfpb > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;	
	%boxcox_survey(
	 data    = _temp2, 
	 subject = subjectid, 
	 var     = pfpb,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x2;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	* tran_lambda=0.1; * = value of macro variable for replicate=0 ;
	tran_paramindex=2;
	varname="pfpb";
	run;
	
	/********************************************************************/
	/*** Dietary constituents #3 - rg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var rg ;
	where rg > 0;
	output out=_min3(keep=minamount) min=minamount ;
	run;
	
	data _min3;
	set _min3;
	minamount=minamount/2;
	tran_paramindex=3;
	length varname $ 32;
	varname="rg";
	run;
	
	data _temp3;
	set inboxcox(keep=subjectid rg &covars.);
	if (rg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp3, 
	 subject = subjectid, 
	 var     = rg,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x3;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.15; * = value of macro variable for replicate=0 ;
	tran_paramindex=3;
	varname="rg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #4 - wg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var wg ;
	where wg > 0;
	output out=_min4(keep=minamount) min=minamount ;
	run;
	
	
	data _min4;
	set _min4;
	minamount=minamount/2;
	tran_paramindex=4;
	length varname $ 32;
	varname="wg";
	run;
	
	data _temp4;
	set inboxcox(keep=subjectid wg &covars.);
	if (wg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp4, 
	 subject = subjectid, 
	 var     = wg,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x4;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.23; * = value of macro variable for replicate=0 ;
	tran_paramindex=4;
	varname="wg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #5 - ssbs: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var ssbs ;
	where ssbs > 0;
	output out=_min5(keep=minamount) min=minamount ;
	run;
	
	
	data _min5;
	set _min5;
	minamount=minamount/2;
	tran_paramindex=5;
	length varname $ 32;
	varname="ssbs";
	run;
	
	data _temp5;
	set inboxcox(keep=subjectid ssbs &covars.);
	if (ssbs > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp5, 
	 subject = subjectid, 
	 var     = ssbs,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x5;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.18;* = value of macro variable for replicate=0 ;
	tran_paramindex=5;
	varname="ssbs";
	run;
	
	/********************************************************************/
	/*** Dietary constituents #6 - vf: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var vf ;
	where vf > 0;
	output out=_min6(keep=minamount) min=minamount ;
	run;
	
	
	data _min6;
	set _min6;
	minamount=minamount/2;
	tran_paramindex=6;
	length varname $ 32;
	varname="vf";
	run;
	
	data _temp6;
	set inboxcox(keep=subjectid vf &covars.);
	if (vf > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp6, 
	 subject = subjectid, 
	 var     = vf,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x6;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=6;
	varname="vf";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #7 - pfab: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var pfab ;
	where pfab > 0;
	output out=_min7(keep=minamount) min=minamount ;
	run;
	
	
	data _min7;
	set _min7;
	minamount=minamount/2;
	tran_paramindex=7;
	length varname $ 32;
	varname="pfab";
	run;
	
	data _temp7;
	set inboxcox(keep=subjectid pfab &covars.);
	if (pfab > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp7, 
	 subject = subjectid, 
	 var     = pfab,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;

	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x7;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.34;* = value of macro variable for replicate=0 ;
	tran_paramindex=7;
	varname="pfab";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #8 - other: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var other ;
	where other > 0;
	output out=_min8(keep=minamount) min=minamount ;
	run;
	
	
	data _min8;
	set _min8;
	minamount=minamount/2;
	tran_paramindex=8;
	length varname $ 32;
	varname="other";
	run;
	
	data _temp8;
	set inboxcox(keep=subjectid other &covars.);
	if (other > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp8, 
	 subject = subjectid, 
	 var     = other,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
		
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x8;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.18;* = value of macro variable for replicate=0 ;
	tran_paramindex=8;
	varname="other";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #9 - water: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var water ;
	where water > 0;
	output out=_min9(keep=minamount) min=minamount ;
	run;
	
	
	data _min9;
	set _min9;
	minamount=minamount/2;
	tran_paramindex=9;
	length varname $ 32;
	varname="water";
	run;
	
	data _temp9;
	set inboxcox(keep=subjectid water &covars.);
	if (water > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp9, 
	 subject = subjectid, 
	 var     = water,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x9;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.41;* = value of macro variable for replicate=0 ;
	tran_paramindex=9;
	varname="water";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #10 - milk: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var milk ;
	where milk > 0;
	output out=_min10(keep=minamount) min=minamount ;
	run;
	
	
	data _min10;
	set _min10;
	minamount=minamount/2;
	tran_paramindex=10;
	length varname $ 32;
	varname="milk";
	run;
	
	data _temp10;
	set inboxcox(keep=subjectid milk &covars.);
	if (milk > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp10, 
	 subject = subjectid, 
	 var     = milk,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x10;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.7;* = value of macro variable for replicate=0 ;
	tran_paramindex=10;
	varname="milk";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #11 - mufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var mufa ;
	where mufa > 0;
	output out=_min11(keep=minamount) min=minamount ;
	run;
	
	data _min11;
	set _min11;
	minamount=minamount/2;
	tran_paramindex=11;
	length varname $ 32;
	varname="mufa";
	run;
	
	data _temp11;
	set inboxcox(keep=subjectid mufa &covars.);
	if (mufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp11, 
	 subject = subjectid, 
	 var     = mufa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x11;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=11;
	varname="mufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #12 - pufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var pufa ;
	where pufa > 0;
	output out=_min12(keep=minamount) min=minamount ;
	run;
	
	data _min12;
	set _min12;
	minamount=minamount/2;
	tran_paramindex=12;
	length varname $ 32;
	varname="pufa";
	run;
	
	data _temp12;
	set inboxcox(keep=subjectid pufa &covars.);
	if (pufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp12, 
	 subject = subjectid, 
	 var     = pufa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x12;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.22;* = value of macro variable for replicate=0 ;
	tran_paramindex=12;
	varname="pufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #13 - sfa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var sfa ;
	where sfa > 0;
	output out=_min13(keep=minamount) min=minamount ;
	run;
	
	data _min13;
	set _min13;
	minamount=minamount/2;
	tran_paramindex=13;
	length varname $ 32;
	varname="sfa";
	run;
	
	data _temp13;
	set inboxcox(keep=subjectid sfa &covars.);
	if (sfa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp13, 
	 subject = subjectid, 
	 var     = sfa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x13;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.29;* = value of macro variable for replicate=0 ;
	tran_paramindex=13;
	varname="sfa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #14 - sugars_free: min. amount and lambda ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var sugars_free ;
	where sugars_free > 0;
	output out=_min14(keep=minamount) min=minamount ;
	run;
	
	
	data _min14;
	set _min14;
	minamount=minamount/2;
	tran_paramindex=14;
	length varname $ 32;
	varname="sugars_free";
	run;
	
	data _temp14;
	set inboxcox(keep=subjectid sugars_free &covars.);
	if (sugars_free > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp14, 
	 subject = subjectid, 
	 var     = sugars_free,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x14;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=14;
	varname="sugars_free";
	run;


	/********************************************************************/
	/*** Dietary constituents #15 - energy: min. amount and lambda    ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum1 min noprint;
	var energy ;
	where energy > 0;
	output out=_min15(keep=minamount) min=minamount ;
	run;
	
	data _min15;
	set _min15;
	minamount=minamount/2;
	tran_paramindex=15;
	length varname $ 32;
	varname="energy";
	run;
	
	data _temp15;
	set inboxcox(keep=subjectid energy &covars.);
	if (energy > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp15, 
	 subject = subjectid, 
	 var     = energy,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x15;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.2;* = value of macro variable for replicate=0 ;
	tran_paramindex=15;
	varname="energy";
	run;

 /*****************************************************************/
 /*** Final step: Append data sets of all dietary constituents  ***/
 /*****************************************************************/
	
	data work.xlambdas_s1_0;
	set _x1-_x15;
	run;
	
	data work.xminamount_s1_0;
	set _min1-_min15;
	label minamount= "Min. non zero amount divided by 2";
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: STANDARDIZE DATA <std_cov_boxcox24hr_conday_minamt>     */
 /*                                                                       */
 /*************************************************************************/

	/* Call the <std_cov_boxcox24hr_conday_minamt> macro with <xlambdas_s1_0> and <xminamount_s1_0> as inputs */
	%std_cov_boxcox24hr_conday_minamt(
	 data                       = _tempstratum1, 
	 prestand_continuous_covars = urine_na_24h , 
	 rec24hr_epis_vars          = plantbev pfpb rg wg ssbs,
	 rec24hr_daily_vars         = vf pfab other milk water mufa pufa sfa sugars_free energy, 
	 boxcox_tran_lambda_data    = xlambdas_s1_0, 
	 minamount_data             = xminamount_s1_0, 
	 print                      = y, 
	 titles                     = 3 ); 
	 
	/* note: contivuous covariates in <prestand_continuous_covars> are standardized and renamed as std_<covar> */

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 1: Fit the multivariate meas. error model (MULTIVAR_MCMC)    */
 /*                                                                       */
 /*************************************************************************/

	/* note: multivar_mcmc may take 48h to 72h to complete */
	
	%multivar_mcmc(
	 data                        = stdcov_stdbc24hr_conday_out, 
	 subject                     = subjectid, 
	 weight_var                  =  , 
	 repeat                      = r24_no, 
	 conday_epis_vars            = conday_plantbev  conday_pfpb  conday_rg  conday_wg  conday_ssbs, 
	 gst_rec24hr_epis_vars       = stdbc_plantbev  stdbc_pfpb  stdbc_rg  stdbc_wg  stdbc_ssbs, 
	 gst_rec24hr_daily_vars      = stdbc_vf  stdbc_pfab  stdbc_other  stdbc_milk  stdbc_water  stdbc_mufa  stdbc_pufa  stdbc_sfa  
								   stdbc_sugars_free  stdbc_energy, 
	 covars_epis_prob            = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 covars_epis_amt             = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 covars_daily_amt            = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 set_seed_mcmc               = 42941, 
	 set_number_mcmc_iterations  = , /*default = 12,000*/
	 set_number_burn_iterations  = , /*default =  2,000*/
	 set_thin                    = , /*default =     25*/
	 prior_sigmau_mean_data      = , 
	 sigmau_constant             = , 
	 gen_inverse                 = y, 
	 print                       = y, 
	 titles                      = 1, 
	 std_print_store             = y, 
	 notes_print                 = y, 
	 out_lib                     = baselib, 
	 out_store_label             = mcmc_s1_rep0, 
	 out_save_label_max5char     = s10, 
	 set_number_saved_out_data   = , 
	 save_mcmc_u_out_data        = y, 
	 set_number_post_mcmc_u_out  = , 
	 traceplots_method1_gpath    = , 
	 traceplots_method2_file_pdf = trace_rep0_s1.pdf, 
	 optional_iml_store_data     = backtran_out, 
	 optional_iml_store_names    = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  
	tran_paramindex tran_lambda tran_center tran_scale minamount 
	 ); 


	/* Save lambdas (Box-Cox transformation lambda values) and minimum amount data */
	data baselib.backtran_out0_s1;
	retain replicate sex;
	set work.backtran_out;
	* indicate replicate number and current stratum value;
	replicate = 0;
	sex = 1;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 2: Simulation of pseudo-individuals (MULTIVAR_DISTRIB)       */
 /*                                                                       */
 /*************************************************************************/

	/* Prepare an input data for the <optional_input_data> option in <multivar_distrib>*/
	proc sort data=_tempstratum1 nodupkey out=optional_input_data(keep= subjectid sex agec_2 
		agec_3 event_primary cens_death urine_na_24h agec);
	by subjectid ;
	run;

	/* Generate pseudo-individuals */
	%multivar_distrib(
	 multivar_mcmc_out_lib           = baselib ,  
	 multivar_mcmc_out_store_label   = mcmc_s1_rep0, 
	 t_weightavg_covariates_list1    = constant1 constant0 constant0 constant0 constant0 constant0 agec_2 agec_3 event_primary 
									   cens_death std_urine_na_24h  , 
	 t_weightavg_covariates_list2    = constant1 constant0 constant0 constant0 constant0 constant1 agec_2 agec_3 event_primary 
									   cens_death std_urine_na_24h  , 
	 set_value_for_weight_cov_list1  = 4, 
	 set_value_for_weight_cov_list2  = 3, 
	 optional_input_data             = optional_input_data , 
	 optional_input_data_var_list    = , 
	 optional_input_mcmc_u_out_data  = , 
	 additional_output_var_list      = sex agec event_primary cens_death urine_na_24h  , 
	 additional_output_subject_var   = subjectid , 
	 output_mcmc_weight_var          = y  , 
	 set_seed_distrib                = 89009890, 
	 set_number_monte_carlo_rand_obs = 500,  
	 print                           = y 
	 ); 

/* Save the Monte Carlo simulation data for current stratum */
	data baselib.mc_t_distrib_out0_s1;
	set mc_t_distrib_out;
	run;

/* delete temporary data sets */
	proc datasets lib=work nolist nodetails;
	delete mc_t_distrib_out optional_input_data stdcov_stdbc24hr_conday_out xlambdas_: 
		   xminamount_: backtran_out ;
	run;
	
	proc datasets lib=baselib nolist nodetails ;
	delete multivar_mcmc_samples_u_outs10 ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: KEEP ONLY PARTICIPANTS IN CURRENT STRATUM (SEX=FEMALES) */
 /*                                                                       */
 /*************************************************************************/

	data _tempstratum2 ;
		set datamrec ;
	* keep only data for current stratum (i.e., females) ;
		if (sex = 2) then output;
	run;
	
	proc sort ;
		by subjectid ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* CURRENT STRATUM = FEMALES ONLY                                        */
 /* PRE NCI STEP: GET MIN AMOUNT AND FIND BEST BOX-COX TRANSFORMATIONS    */
 /*                                                                       */
 /*************************************************************************/

	/* Output data from the first 24-h dietary recall completed only */
	proc sort data = _tempstratum2 nodupkey out=inboxcox;
	by subjectid ;
	run;
	
	/* make a list of covariates for consistency */
	%let covars  = seq2 seq3 seq4 seq5 r24_wkd urine_na_24h agec_2 agec_3 event_primary cens_death ;
	
	/* ensure that the macro variable <best_lambda> is available outside the <boxcox_survey> macro */
	%global best_lambda ;

	/********************************************************************/
	/*** Dietary constituents #1 - plantbev: min. amount and lambda   ***/
	/********************************************************************/

	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var plantbev ;
	where plantbev > 0;
	output out=_min1(keep=minamount) min=minamount ;
	run;
	
	data _min1;
	set _min1;
	minamount=minamount/2;
	tran_paramindex=1;
	length varname $ 32;
	varname="plantbev";
	run;

	data _temp1;
	set inboxcox(keep=subjectid plantbev &covars);
	if (plantbev > 0) then output;
	run;

	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp1, 
	 subject = subjectid, 
	 var     = plantbev,  
	 covars  = &covars. ,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;

	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x1;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	* tran_lambda=0.29; * = value of macro variable for replicate=0, stratum=females ;
	tran_paramindex=1;
	varname="plantbev";
	run;

	/* note: steps 1-3 are repeated for all other dietary constituents of the HEFI-2019 below */

	/********************************************************************/
	/*** Dietary constituents #2 - pfpb: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var pfpb ;
	where pfpb > 0;
	output out=_min2(keep=minamount) min=minamount ;
	run;
	
	data _min2;
	set _min2;
	minamount=minamount/2;
	tran_paramindex=2;
	length varname $ 32;
	varname="pfpb";
	run;
	
	data _temp2;
	set inboxcox(keep=subjectid pfpb &covars. );
	if (pfpb > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;	
	%boxcox_survey(
	 data    = _temp2, 
	 subject = subjectid, 
	 var     = pfpb,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x2;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	* tran_lambda=0.1; * = value of macro variable for replicate=0 ;
	tran_paramindex=2;
	varname="pfpb";
	run;
	
	/********************************************************************/
	/*** Dietary constituents #3 - rg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var rg ;
	where rg > 0;
	output out=_min3(keep=minamount) min=minamount ;
	run;
	
	data _min3;
	set _min3;
	minamount=minamount/2;
	tran_paramindex=3;
	length varname $ 32;
	varname="rg";
	run;
	
	data _temp3;
	set inboxcox(keep=subjectid rg &covars.);
	if (rg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp3, 
	 subject = subjectid, 
	 var     = rg,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x3;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.15; * = value of macro variable for replicate=0 ;
	tran_paramindex=3;
	varname="rg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #4 - wg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var wg ;
	where wg > 0;
	output out=_min4(keep=minamount) min=minamount ;
	run;
	
	
	data _min4;
	set _min4;
	minamount=minamount/2;
	tran_paramindex=4;
	length varname $ 32;
	varname="wg";
	run;
	
	data _temp4;
	set inboxcox(keep=subjectid wg &covars.);
	if (wg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp4, 
	 subject = subjectid, 
	 var     = wg,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x4;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.23; * = value of macro variable for replicate=0 ;
	tran_paramindex=4;
	varname="wg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #5 - ssbs: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var ssbs ;
	where ssbs > 0;
	output out=_min5(keep=minamount) min=minamount ;
	run;
	
	
	data _min5;
	set _min5;
	minamount=minamount/2;
	tran_paramindex=5;
	length varname $ 32;
	varname="ssbs";
	run;
	
	data _temp5;
	set inboxcox(keep=subjectid ssbs &covars.);
	if (ssbs > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp5, 
	 subject = subjectid, 
	 var     = ssbs,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x5;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.18;* = value of macro variable for replicate=0 ;
	tran_paramindex=5;
	varname="ssbs";
	run;
	
	/********************************************************************/
	/*** Dietary constituents #6 - vf: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var vf ;
	where vf > 0;
	output out=_min6(keep=minamount) min=minamount ;
	run;
	
	
	data _min6;
	set _min6;
	minamount=minamount/2;
	tran_paramindex=6;
	length varname $ 32;
	varname="vf";
	run;
	
	data _temp6;
	set inboxcox(keep=subjectid vf &covars.);
	if (vf > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp6, 
	 subject = subjectid, 
	 var     = vf,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x6;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=6;
	varname="vf";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #7 - pfab: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var pfab ;
	where pfab > 0;
	output out=_min7(keep=minamount) min=minamount ;
	run;
	
	
	data _min7;
	set _min7;
	minamount=minamount/2;
	tran_paramindex=7;
	length varname $ 32;
	varname="pfab";
	run;
	
	data _temp7;
	set inboxcox(keep=subjectid pfab &covars.);
	if (pfab > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp7, 
	 subject = subjectid, 
	 var     = pfab,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;

	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x7;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.34;* = value of macro variable for replicate=0 ;
	tran_paramindex=7;
	varname="pfab";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #8 - other: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var other ;
	where other > 0;
	output out=_min8(keep=minamount) min=minamount ;
	run;
	
	
	data _min8;
	set _min8;
	minamount=minamount/2;
	tran_paramindex=8;
	length varname $ 32;
	varname="other";
	run;
	
	data _temp8;
	set inboxcox(keep=subjectid other &covars.);
	if (other > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp8, 
	 subject = subjectid, 
	 var     = other,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
		
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x8;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.18;* = value of macro variable for replicate=0 ;
	tran_paramindex=8;
	varname="other";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #9 - water: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var water ;
	where water > 0;
	output out=_min9(keep=minamount) min=minamount ;
	run;
	
	
	data _min9;
	set _min9;
	minamount=minamount/2;
	tran_paramindex=9;
	length varname $ 32;
	varname="water";
	run;
	
	data _temp9;
	set inboxcox(keep=subjectid water &covars.);
	if (water > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp9, 
	 subject = subjectid, 
	 var     = water,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x9;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.41;* = value of macro variable for replicate=0 ;
	tran_paramindex=9;
	varname="water";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #10 - milk: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var milk ;
	where milk > 0;
	output out=_min10(keep=minamount) min=minamount ;
	run;
	
	
	data _min10;
	set _min10;
	minamount=minamount/2;
	tran_paramindex=10;
	length varname $ 32;
	varname="milk";
	run;
	
	data _temp10;
	set inboxcox(keep=subjectid milk &covars.);
	if (milk > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp10, 
	 subject = subjectid, 
	 var     = milk,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x10;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.7;* = value of macro variable for replicate=0 ;
	tran_paramindex=10;
	varname="milk";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #11 - mufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var mufa ;
	where mufa > 0;
	output out=_min11(keep=minamount) min=minamount ;
	run;
	
	data _min11;
	set _min11;
	minamount=minamount/2;
	tran_paramindex=11;
	length varname $ 32;
	varname="mufa";
	run;
	
	data _temp11;
	set inboxcox(keep=subjectid mufa &covars.);
	if (mufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp11, 
	 subject = subjectid, 
	 var     = mufa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x11;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=11;
	varname="mufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #12 - pufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var pufa ;
	where pufa > 0;
	output out=_min12(keep=minamount) min=minamount ;
	run;
	
	data _min12;
	set _min12;
	minamount=minamount/2;
	tran_paramindex=12;
	length varname $ 32;
	varname="pufa";
	run;
	
	data _temp12;
	set inboxcox(keep=subjectid pufa &covars.);
	if (pufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp12, 
	 subject = subjectid, 
	 var     = pufa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x12;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.22;* = value of macro variable for replicate=0 ;
	tran_paramindex=12;
	varname="pufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #13 - sfa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var sfa ;
	where sfa > 0;
	output out=_min13(keep=minamount) min=minamount ;
	run;
	
	data _min13;
	set _min13;
	minamount=minamount/2;
	tran_paramindex=13;
	length varname $ 32;
	varname="sfa";
	run;
	
	data _temp13;
	set inboxcox(keep=subjectid sfa &covars.);
	if (sfa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp13, 
	 subject = subjectid, 
	 var     = sfa,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x13;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.29;* = value of macro variable for replicate=0 ;
	tran_paramindex=13;
	varname="sfa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #14 - sugars_free: min. amount and lambda ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var sugars_free ;
	where sugars_free > 0;
	output out=_min14(keep=minamount) min=minamount ;
	run;
	
	
	data _min14;
	set _min14;
	minamount=minamount/2;
	tran_paramindex=14;
	length varname $ 32;
	varname="sugars_free";
	run;
	
	data _temp14;
	set inboxcox(keep=subjectid sugars_free &covars.);
	if (sugars_free > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp14, 
	 subject = subjectid, 
	 var     = sugars_free,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x14;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=14;
	varname="sugars_free";
	run;


	/********************************************************************/
	/*** Dietary constituents #15 - energy: min. amount and lambda    ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=_tempstratum2 min noprint;
	var energy ;
	where energy > 0;
	output out=_min15(keep=minamount) min=minamount ;
	run;
	
	data _min15;
	set _min15;
	minamount=minamount/2;
	tran_paramindex=15;
	length varname $ 32;
	varname="energy";
	run;
	
	data _temp15;
	set inboxcox(keep=subjectid energy &covars.);
	if (energy > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp15, 
	 subject = subjectid, 
	 var     = energy,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 );
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */	
	data _x15;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	*tran_lambda=0.2;* = value of macro variable for replicate=0 ;
	tran_paramindex=15;
	varname="energy";
	run;

 /*****************************************************************/
 /*** Final step: Append data sets of all dietary constituents  ***/
 /*****************************************************************/
	
	data work.xlambdas_s2_0;
	set _x1-_x15;
	run;
	
	data work.xminamount_s2_0;
	set _min1-_min15;
	label minamount= "Min. non zero amount divided by 2";
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: STANDARDIZE DATA <std_cov_boxcox24hr_conday_minamt>     */
 /*                                                                       */
 /*************************************************************************/

	/* Call the <std_cov_boxcox24hr_conday_minamt> macro with <xlambdas_s2_0> and <xminamount_s2_0> as inputs */
	%std_cov_boxcox24hr_conday_minamt(
	 data                       = _tempstratum2, 
	 prestand_continuous_covars = urine_na_24h , 
	 rec24hr_epis_vars          = plantbev pfpb rg wg ssbs,
	 rec24hr_daily_vars         = vf pfab other milk water mufa pufa sfa sugars_free energy, 
	 boxcox_tran_lambda_data    = xlambdas_s2_0, 
	 minamount_data             = xminamount_s2_0, 
	 print                      = y, 
	 titles                     = 3 ); 
	 
	/* note: contivuous covariates in <prestand_continuous_covars> are standardized and renamed as std_<covar> */

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 1: Fit the multivariate meas. error model (MULTIVAR_MCMC)    */
 /*                                                                       */
 /*************************************************************************/

	/* note: multivar_mcmc may take 48h to 72h to complete */
	
	%multivar_mcmc(
	 data                        = stdcov_stdbc24hr_conday_out, 
	 subject                     = subjectid, 
	 weight_var                  =  , 
	 repeat                      = r24_no, 
	 conday_epis_vars            = conday_plantbev  conday_pfpb  conday_rg  conday_wg  conday_ssbs, 
	 gst_rec24hr_epis_vars       = stdbc_plantbev  stdbc_pfpb  stdbc_rg  stdbc_wg  stdbc_ssbs, 
	 gst_rec24hr_daily_vars      = stdbc_vf  stdbc_pfab  stdbc_other  stdbc_milk  stdbc_water  stdbc_mufa  stdbc_pufa  stdbc_sfa  
								   stdbc_sugars_free  stdbc_energy, 
	 covars_epis_prob            = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 covars_epis_amt             = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 covars_daily_amt            = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  , 
	 set_seed_mcmc               = 42941, 
	 set_number_mcmc_iterations  = , /*default = 12,000*/
	 set_number_burn_iterations  = , /*default =  2,000*/
	 set_thin                    = , /*default =     25*/
	 prior_sigmau_mean_data      = , 
	 sigmau_constant             = , 
	 gen_inverse                 = y, 
	 print                       = y, 
	 titles                      = 1, 
	 std_print_store             = y, 
	 notes_print                 = y, 
	 out_lib                     = baselib, 
	 out_store_label             = mcmc_s2_rep0, 
	 out_save_label_max5char     = s20, 
	 set_number_saved_out_data   = , 
	 save_mcmc_u_out_data        = y, 
	 set_number_post_mcmc_u_out  = , 
	 traceplots_method1_gpath    = , 
	 traceplots_method2_file_pdf = trace_rep0_s2.pdf, 
	 optional_iml_store_data     = backtran_out, 
	 optional_iml_store_names    = constant1 seq2 seq3 seq4 seq5 r24_wkd agec_2 agec_3 event_primary cens_death std_urine_na_24h  
	tran_paramindex tran_lambda tran_center tran_scale minamount 
	 ); 


	/* Save lambdas (Box-Cox transformation lambda values) and minimum amount data */
	data baselib.backtran_out0_s2;
	retain replicate sex;
	set work.backtran_out;
	* indicate replicate number and current stratum value;
	replicate = 0;
	sex = 2;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 2: Simulation of pseudo-individuals (MULTIVAR_DISTRIB)       */
 /*                                                                       */
 /*************************************************************************/

	/* Prepare an input data for the <optional_input_data> option in <multivar_distrib>*/
	proc sort data=_tempstratum2 nodupkey out=optional_input_data(keep= subjectid sex agec_2 
		agec_3 event_primary cens_death urine_na_24h agec);
	by subjectid ;
	run;

	/* Generate pseudo-individuals */
	%multivar_distrib(
	 multivar_mcmc_out_lib           = baselib ,  
	 multivar_mcmc_out_store_label   = mcmc_s2_rep0, 
	 t_weightavg_covariates_list1    = constant1 constant0 constant0 constant0 constant0 constant0 agec_2 agec_3 event_primary 
									   cens_death std_urine_na_24h  , 
	 t_weightavg_covariates_list2    = constant1 constant0 constant0 constant0 constant0 constant1 agec_2 agec_3 event_primary 
									   cens_death std_urine_na_24h  , 
	 set_value_for_weight_cov_list1  = 4, 
	 set_value_for_weight_cov_list2  = 3, 
	 optional_input_data             = optional_input_data , 
	 optional_input_data_var_list    = , 
	 optional_input_mcmc_u_out_data  = , 
	 additional_output_var_list      = sex agec event_primary cens_death urine_na_24h  , 
	 additional_output_subject_var   = subjectid , 
	 output_mcmc_weight_var          = y  , 
	 set_seed_distrib                = 89009890, 
	 set_number_monte_carlo_rand_obs = 500,  
	 print                           = y 
	 ); 

/* Save the Monte Carlo simulation data for current stratum */
	data baselib.mc_t_distrib_out0_s2;
	set mc_t_distrib_out;
	run;

/* delete temporary data sets */
	proc datasets lib=work nolist nodetails;
	delete mc_t_distrib_out optional_input_data stdcov_stdbc24hr_conday_out xlambdas_: 
		   xminamount_: backtran_out ;
	run;
	
	proc datasets lib=baselib nolist nodetails ;
	delete multivar_mcmc_samples_u_outs20 ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: PREPARE MC_T_DISTRIB_OUT DATA FOR FURTHER ANALYSIS     */
 /*                                                                       */
 /*************************************************************************/

	data baselib.usintake_mc_t0;
		set baselib.mc_t_distrib_out0_s1-baselib.mc_t_distrib_out0_s2;
	
		/* note: Monte Carlo simulations from both stratum are appended */
	
	* Divide weights by <set_number_monte_carlo_rand_obs> value (i.e., number of pseudo-individuals);

		weight_nw_sumw_div = weight_nw_sumw / 500;
	
	* Rename dietary constituents ;
	
		array outmc (*) mc_t1-mc_t15 ;
		array clean (*) plantbev pfpb rg wg ssbs vf pfab other milk water mufa pufa sfa sugars_free energy ;
		do i=1 to dim(outmc);
			clean(i)=outmc(i);
		end;
	
	run;

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: Apply the HEFI-2019 scoring algorithm                  */
 /*                                                                       */
 /*************************************************************************/

 /* Call the <HEFI2019> macro to apply the HEFI-2019 scoring algorithm*/
	%HEFI2019(
	 indata             = baselib.usintake_mc_t0, 
	 vegfruits          = vf, 
	 wholegrfoods       = wg, 
	 nonwholegrfoods    = rg, 
	 profoodsanimal     = pfab, 
	 profoodsplant      = pfpb, 
	 otherfoods         = other, 
	 waterhealthybev    = water, 
	 unsweetmilk        = milk, 
	 unsweetplantbevpro = plantbev, 
	 otherbeverages     = ssbs, 
	 mufat              = mufa, 
	 pufat              = pufa, 
	 satfat             = sfa, 
	 freesugars         = sugars_free, 
	 sodium             = urine_na_24h, 
	 energy             = energy, 
	 outdata            = baselib.usintake_mc_t0
	 ); 

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: OUTPUT USUAL INTAKE DISTRIBUTION BASED SIMULATIONS     */
 /*                                                                       */
 /*************************************************************************/

	/* note: distribution are generated using the NCI macro <percentiles_Survey> */

/* Total HEFI-2019 score - full eligible sample */
    %percentiles_Survey(data      = baselib.usintake_mc_t0,
                        byvar     = ,
                        var       = HEFI2019_TOTAL_SCORE ,
                        weight    = weight_nw_sumw_div,
                        cutpoints = ,
                        print     = N,
                        ntitle    = 3
                        );
	
	/* Format and save output data */
		data reslib.distribtotal_w0 ;
		retain replicate varname;
		set _percentiles;
		* add variable identifier ;
			replicate=0;
			length varname $ 32;
			varname="HEFI2019_TOTAL_SCORE";
		run;

/* Total HEFI-2019 score - by censoring indicator */
    %percentiles_Survey(data      = baselib.usintake_mc_t0,
                        byvar     = cens_death,
                        var       = HEFI2019_TOTAL_SCORE ,
                        weight    = weight_nw_sumw_div,
                        cutpoints = ,
                        print     = N,
                        ntitle    = 3
                        );
	
	/* Format and save output data */
		data reslib.distribtotal_w_cens_death0 ;
		retain replicate varname cens_death;
		set _percentiles;
		* add variable identifier ;
			replicate=0;
			length varname $ 32;
			varname="HEFI2019_TOTAL_SCORE";
		run;

/* Energy intake - full eligible sample */
    %percentiles_Survey(data      = baselib.usintake_mc_t0,
                        byvar     = ,
                        var       = energy  ,
                        weight    = weight_nw_sumw_div,
                        cutpoints = ,
                        print     = N,
                        ntitle    = 3
                        );
	
	/* Format and save output data */
		data reslib.distribraw_w0 ;
		retain replicate varname;
		set _percentiles;
		* add variable identifier ;
			replicate=0;
			length varname $ 32;
			varname="energy";
		run;

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: ADDITIONAL DATA WITH PERCENTILE KNOTS VALUE FOR SPLINE */
 /*                                                                       */
 /*************************************************************************/

	data _total;
	retain food ;
		set reslib.distribtotal_w0(keep=varname /* spline 4-knot: */ Pctile5 Pctile35 Pctile65 Pctile95 rename=(varname=food));
	rename  pctile5=k1 pctile35=k2 pctile65=k3 pctile95=k4 ;
	run;
	
	data _food;
		set reslib.distribraw_w0(keep=varname /* spline 4-knot: */ Pctile5 Pctile35 Pctile65 Pctile95 rename=(varname=food));
	rename pctile5=k1 pctile35=k2 pctile65=k3 pctile95=k4  ;
	run;
	
	data nci.hefi_spl;
		set _total _food;
	run;
	

 /*************************************************************************/
 /*                                                                       */
 /*                              END OF CODE                              */
 /*                                                                       */
 /*************************************************************************/
