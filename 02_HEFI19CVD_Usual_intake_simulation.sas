 /*************************************************************************/
 /*                                                                       */
 /*               HEFI-2019 and CVD risk in the UK Biobank                */
 /*                                                                       */
 /*                  Brassard et al. Am J Clin Nutr 2022                  */
 /*       Code 2: Usual intake simulation at the participant level        */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                              Version 1.1                              */
 /*                               25AUG2022                               */
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
	%let suffix=_f ;

/* National Cancer Institute (NCI) macros */
	/* boxcox svy */ %include "&path./Macros/boxcox_survey.macro.v1.2.sas";
	/* std boxcox */ %include "&path./Macros/std_cov_boxcox24hr_conday_minamt_macro_v2.0.sas";
	/* multi MCMC */ %include "&path./Macros/multivar_mcmc_macro_v2.1.sas";
	/* multi dist */ %include "&path./Macros/multivar_distrib_macro_v2.1.sas";
	
	/* Available at: https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0 */
	
/* HEFI-2019 scoring algorithm macro */
	/* hefi-2019 */ %include "&path./Macros/hefi2019.scoring.macro.sas";
	
	/* Available at: https://github.com/didierbrassard/hefi2019/ */
	
/* create or indicate library folder */
	options dlcreatedir;
	libname temp "&path./Temp/";
	libname fmtdata "&path./Fmtdata/";
	libname distrib "&path./NCI/MCMC_distrib/Results/";  /* data of population-level distribution of usual intakes (not shown) */
	
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

/* list of covariates */
%let covar= ageat24hr deprivation bmi0 sedentary0 met_level0 
	 sexmeno university employ0 region british famhist0 
	 roh0 smk0 physact0 dietchange0 dietsuppl rx0 hrt0 riskfactor ;

	data _preNCI(drop=valid24hr included rename=(_roh0=roh0 nomat=subjectid));
		merge fmtdata.diet24hr(in=foods drop=mean:)
			  fmtdata.valid24hr(in=valid keep=nomat valid24hr where=(valid24hr=1))
			  fmtdata.date24hr( keep=nomat r24_date rename=(r24_date=date_24hr))
			  fmtdata.screening(keep=nomat sex included cens_death event_date: cens_primary event_primary time_to_primary )
			  fmtdata.urine_na_24h(keep=nomat urine_na_24h)
			  fmtdata.postimp(keep=nomat _imputation_ &covar missing_top2 where=(_Imputation_=1))
			  ;
		by nomat ;
		
	**************************************************;
	* REGION: collapse similar  regions for analysis  ;
	**************************************************;
	
	length region2 $ 8;
		if 1<=region<=8 then region2='England';
		else if region=9 then region2='Wales';
		else if region=10 then region2='Scotland';

	**************************************************;
	* PLANTBEV: distribute intake to beverage / pfpb  ;
	**************************************************;

	* Soy beverages are distributed to beverages and plant-based protein foods;
	array pfpb(5) pfpb1-pfpb5;
	array water(5) water1-water5;
	array plantbev(5) plantbev1-plantbev5;
	
	do i=1 to dim(plantbev);
		* add ml of plantbev to water (and other healthy beverage);
		if not missing(plantbev(i)) then water(i) = water(i) + plantbev(i);
		
		* add RA of plantbev to pfpb;
		if not missing(plantbev(i)) then pfpb(i) = (plantbev(i)/250) + pfpb(i);
	end;
	
	drop i;

	**************************************************;
	* ROH USE: collapse uncommon levels for analysis  ;
	**************************************************;

	if roh0 in (0 1) then _roh0 = 0;* approx. 6 percent;
		else if roh0 = 2 then _roh0 = 1; * approx. 94 percent;
	
	drop roh0;
	
	label _roh0 = "Drinker indicator (f.20117)";
	
	*******************************************************;
	* Output valid 24hr, eligible participants, no missing ;
	*******************************************************;
		
	if valid & included=1 & missing_top2=0 then output;
	run;
	
		/* note: eligible sample = 136,698 obs */
	
	proc sort data=_preNCI;
		by subjectid ;
	run;
		
	/* Use <proc glmmod> to output dummy-coded variables for categorical covariates */
	proc glmmod data=_preNCI outdesign=_GLMDesign outparm=_GLMParm noprint;
		class sexmeno region2 famhist0 employ0 university smk0 roh0 rx0 hrt0 dietsuppl 
			riskfactor physact0 sex;
		model subjectid = sexmeno region2 famhist0 employ0 university smk0 roh0 rx0 hrt0 
			dietsuppl riskfactor physact0 sex;
	run;
	
	
	/* output list of level names for dummy-coded variables based on <_GLMParm> */
	data _GLMParm(keep=EffName Level);
		set _GLMParm(where=(EffName not in ("Intercept" "Constante")));
	* some formatting needed for efficient renaming ;
		Length Level $ 32;
		Level=cats(sexmeno, region2, famhist0, employ0, university, smk0, roh0, rx0, 
			hrt0, dietsuppl, riskfactor, physact0, sex);
		Level=Compress(Level, ".");
		Level=cats(EffName, "_", Level);
	run;
	
	Proc Sql noprint;
		Select Level into : Lvl_List separated by ' ' from _GLMParm;
	Quit;
	
	/* Rename dummy-coded variables according to ordered levels */
	data _dummyz (drop=Col2-Col36 ith);
		set _GLMDesign(Keep=subjectid Col2-Col36);
		array old(*) Col2-Col36;
		array new(*) &Lvl_list;
		
		* note: <Lvl_list> = sexmeno_1 sexmeno_2 sexmeno_3 sexmeno_4 region2_England region2_Scotland region2_Wales 
	 famhist0_0 famhist0_1 famhist0_2 employ0_1 employ0_2 employ0_3 university_0 university_1 smk0_0 smk0_1 smk0_2 roh0_0 roh0_1 rx0_0 
	 rx0_1 rx0_2 hrt0_0 hrt0_1 dietsuppl_0 dietsuppl_1 riskfactor_0 riskfactor_1 riskfactor_2 physact0_1 physact0_2 physact0_3 
	 sex_Female sex_Male ;
	
		do ith=1 to dim(old);
			new(ith)=old(ith);
		end;
	run;

	/* combine original data with dummy-coded covariates */
	data fmtdata.preNCI;
		merge _preNCI(in=a) _dummyz;
			by subjectid;
	* output participants in preNCI;
		if a then output;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: CREATE BOOTSTRAP SAMPLES                                */
 /*                                                                       */
 /*************************************************************************/

	data _original;
	retain replicate;
	set fmtdata.preNCI;
	
	* make replicate data id for original sample;
	replicate=0;
	
	keep /* identification: */
		subjectid replicate
		
		/* 24-h dietary recall specific variables: */
		vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 sugars_free1 energy1 pfpb1 rg1 wg1 ssbs1
		vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 sugars_free2 energy2 pfpb2 rg2 wg2 ssbs2
		vf3 pfab3 other3 water3 milk3 mufa3 pufa3 sfa3 sugars_free3 energy3 pfpb3 rg3 wg3 ssbs3
		vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 sugars_free4 energy4 pfpb4 rg4 wg4 ssbs4 
		vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 sugars_free5 energy5 pfpb5 rg5 wg5 ssbs5
		
		/* Covariates, categorical (reference level for each are not kept): */
		sexmeno_2 sexmeno_3 sexmeno_4 region2_Scotland region2_Wales famhist0_1 famhist0_2 employ0_2 employ0_3 university_1
		smk0_1 smk0_2 roh0_1 rx0_1 rx0_2 hrt0_1 dietsuppl_1 riskfactor_1 riskfactor_2 physact0_2 physact0_3
		
		/* Covariates, continuous: */
		ageat24hr deprivation sedentary0 bmi0 urine_na_24h
		;
	run;
	
	/* perfom simple random sampling with replacement 250 times */
	proc surveyselect data=_original(drop=replicate) out=_workbootrep seed=1563783 
			method=URS sampsize=136698 outhits reps=250;
	run;
	 
	/* append original data with bootstrap data */
	data baselib.original_and_boot;
	set _original(in=a) _workbootrep ;
	if a then Replicate=0;
	run;
	 
	proc sort data=baselib.original_and_boot;
	by replicate;
	run;
	
	/* clean temporary data */
	proc datasets lib=work nolist nodetails;
	delete _original _workbootrep;
	run;

	/* create new subject id variable, by bootstrap sample */
	data baselib.original_and_boot;
	set baselib.original_and_boot;
	by replicate;
	retain replicaterowid;
	if first.replicate then replicaterowid=1;
	else replicaterowid=replicaterowid + 1;
	run;

	data baselib.original_and_boot(drop=numberhits );
	retain replicate replicaterowid subjectid ;
	set baselib.original_and_boot;
	run;
	 
	proc sort data=baselib.original_and_boot;
	by replicate replicaterowid ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: PREPARE DATA SETS NEEDED FOR MULTIVARIATE METHOD        */
 /*                                                                       */
 /*************************************************************************/

	
	 /*************************************************************************/
	 /*                                                                       */
	 /*                         BOOTSTRAP LOOP START                          */
	 /*     THE FOLLOWING ANALYSES ARE REQUIRED FOR EACH BOOTSTRAP SAMPLE     */
	 /*                                                                       */
	 /*************************************************************************/

	%macro replicateloop(replicfirst=,repliclast=) ;
	options cpucount=actual;
	
	/* Begin the replicate loop */
	%do replicnum=&replicfirst %to &repliclast;

	/* Transfer SAS log to file */
	filename UILOG "&path./NCI/MCMC&suffix/nci_epi_rep&replicfirst._to_&repliclast..log";
	proc printto log = UILOG; run;
	
	/* Print contextual information to log */
	%put # Started on: &sysdate9 (&systime) ;
	%put # Current analysis is for &=replicnum / &repliclast;
	
	/* Define output library according to value of replicnum */
		%if (&replicnum=0) %then %let outlib=baselib;
		%else %let outlib=bootlib;
	
	/************************************************/
	/*    Begin analysis for current <replicnum>    */
	/************************************************/

	/* Make data with 1-row per subject (needed AFTER NCI) for current replicate, based on original_and_boot data */
	data subj1rec;
	set baselib.original_and_boot (where = (replicate = &replicnum.));
	* constant1 - used to have an intercept in NCI multivariate models;
	constant1=1;
	run;

	proc sort data=subj1rec;
	by replicaterowid;
	run;

	/* Make data with m-rows per subject where m = number of 24-h recall (needed FOR NCI)*/
	data datamrec (drop=vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 sugars_free1 
			energy1 pfpb1 rg1 wg1 ssbs1 vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 
			sugars_free2 energy2 pfpb2 rg2 wg2 ssbs2 vf3 pfab3 other3 water3 milk3 mufa3 
			pufa3 sfa3 sugars_free3 energy3 pfpb3 rg3 wg3 ssbs3 vf4 pfab4 other4 water4 
			milk4 mufa4 pufa4 sfa4 sugars_free4 energy4 pfpb4 rg4 wg4 ssbs4 vf5 pfab5 
			other5 water5 milk5 mufa5 pufa5 sfa5 sugars_free5 energy5 pfpb5 rg5 wg5 ssbs5 
			ith);
		* drop recall-specific names ;
		retain replicate replicaterowid r24_no;
		set subj1rec;
		
		* generic variable names ;
		array foodlist0 (*) vf pfab other water milk mufa pufa sfa sugars_free energy 
			pfpb rg wg ssbs;
			
		* assign generic variable names to recall specific variable names ;
		
		* Recall #1: output data if there are no missing values  ;
		if nmiss(of vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 sugars_free1 
			energy1 pfpb1 rg1 wg1 ssbs1)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr1 (*) vf1 pfab1 other1 water1 milk1 mufa1 pufa1 sfa1 
					sugars_free1 energy1 pfpb1 rg1 wg1 ssbs1;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr1(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=1;
				* writes observation for current recall ;
				output;
			end;

		* Recall #2: output data if there are no missing values  ;
		if nmiss(of vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 sugars_free2 
			energy2 pfpb2 rg2 wg2 ssbs2)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr2 (*) vf2 pfab2 other2 water2 milk2 mufa2 pufa2 sfa2 
					sugars_free2 energy2 pfpb2 rg2 wg2 ssbs2;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr2(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=2;
				* writes observation for current recall ;
				output;
			end;
			
		* Recall #3: output data if there are no missing values  ;
		if nmiss(of vf3 pfab3 other3 water3 milk3 mufa3 pufa3 sfa3 sugars_free3 
			energy3 pfpb3 rg3 wg3 ssbs3)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr3 (*) vf3 pfab3 other3 water3 milk3 mufa3 pufa3 sfa3 
					sugars_free3 energy3 pfpb3 rg3 wg3 ssbs3;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr3(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=3;
				* writes observation for current recall ;
				output;
			end;

		* Recall #4: output data if there are no missing values  ;
		if nmiss(of vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 sugars_free4 
			energy4 pfpb4 rg4 wg4 ssbs4)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr4 (*) vf4 pfab4 other4 water4 milk4 mufa4 pufa4 sfa4 
					sugars_free4 energy4 pfpb4 rg4 wg4 ssbs4;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr4(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=4;
				* writes observation for current recall ;
				output;
			end;
			
		* Recall #5: output data if there are no missing values  ;
		if nmiss(of vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 sugars_free5 
			energy5 pfpb5 rg5 wg5 ssbs5)=0 then do;
				* change recall-specific to generic variable names;
				array foodlistr5 (*) vf5 pfab5 other5 water5 milk5 mufa5 pufa5 sfa5 
					sugars_free5 energy5 pfpb5 rg5 wg5 ssbs5;
	
				do ith=1 to dim(foodlist0);
					foodlist0(ith)=foodlistr5(ith);
					* reads as: for the <ith> dietary constituent, generic name = recall-specific ... ;
				end;
				* sequential recall id ;
				r24_no=5;
				* writes observation for current recall ;
				output;
			end;
	run;

	proc sort data=datamrec;
		by replicaterowid ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: GET MIN AMOUNT AND FIND BEST BOX-COX TRANSFORMATIONS    */
 /*                                                                       */
 /*************************************************************************/

	/* Output data from the first 24-h dietary recall completed only */
	proc sort data = datamrec nodupkey out=inboxcox;
		by replicaterowid ;
	run;
	
	/* make a list of covariates for consistency */
	%let covars  = ageat24hr deprivation sedentary0 bmi0 urine_na_24h sexmeno_2 
		sexmeno_3 sexmeno_4 region2_Scotland region2_Wales famhist0_1 famhist0_2 employ0_2 employ0_3 university_1
		smk0_1 smk0_2 roh0_1 rx0_1 rx0_2 hrt0_1 dietsuppl_1 riskfactor_1 riskfactor_2 physact0_2 physact0_3 ;
	
	/* ensure that the macro variable <best_lambda> is available outside the <boxcox_survey> macro */
	%global best_lambda ;

	/********************************************************************/
	/*** Dietary constituents #1 - pfpb: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var pfpb ;
	where pfpb > 0;
	output out=_min1(keep=minamount) min=minamount ;
	run;
	
	data _min1;
	set _min1;
	minamount=minamount/2;
	tran_paramindex=1;
	length varname $ 32;
	varname="pfpb";
	run;
	
	data _temp1;
	set inboxcox(keep=replicaterowid pfpb &covars. );
	if (pfpb > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;	
	%boxcox_survey(
	 data    = _temp1, 
	 subject = replicaterowid, 
	 var     = pfpb,  
	 covars  = &covars.,  
	 weight  = ,  
	 print   = N, 
	 plot    = N, 
	 ntitle  = 4 ); 
	options notes source;
	
	/* 3) Save <best_lambda> macro variable from boxcox_survey for current dietary constituent */
	data _x1;
	length varname $ 32;
	tran_lambda=&best_lambda ; * = macro variable reference ;
	* tran_lambda=0.04; * = value of macro variable for replicate=0 ;
	tran_paramindex=1;
	varname="pfpb";
	run;
	
	/* note: steps 1-3 are repeated for all other dietary constituents of the HEFI-2019 below */
	
	/********************************************************************/
	/*** Dietary constituents #2 - rg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var rg ;
	where rg > 0;
	output out=_min2(keep=minamount) min=minamount ;
	run;
	
	data _min2;
	set _min2;
	minamount=minamount/2;
	tran_paramindex=2;
	length varname $ 32;
	varname="rg";
	run;
	
	data _temp2;
	set inboxcox(keep=replicaterowid rg &covars.);
	if (rg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp2, 
	 subject = replicaterowid, 
	 var     = rg,  
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
	*tran_lambda=0.09; * = value of macro variable for replicate=0 ;
	tran_paramindex=2;
	varname="rg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #3 - wg: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var wg ;
	where wg > 0;
	output out=_min3(keep=minamount) min=minamount ;
	run;
	
	
	data _min3;
	set _min3;
	minamount=minamount/2;
	tran_paramindex=3;
	length varname $ 32;
	varname="wg";
	run;
	
	data _temp3;
	set inboxcox(keep=replicaterowid wg &covars.);
	if (wg > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp3, 
	 subject = replicaterowid, 
	 var     = wg,  
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
	*tran_lambda=0.09; * = value of macro variable for replicate=0 ;
	tran_paramindex=3;
	varname="wg";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #4 - ssbs: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var ssbs ;
	where ssbs > 0;
	output out=_min4(keep=minamount) min=minamount ;
	run;
	
	
	data _min4;
	set _min4;
	minamount=minamount/2;
	tran_paramindex=4;
	length varname $ 32;
	varname="ssbs";
	run;
	
	data _temp4;
	set inboxcox(keep=replicaterowid ssbs &covars.);
	if (ssbs > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp4, 
	 subject = replicaterowid, 
	 var     = ssbs,  
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
	*tran_lambda=0.22;* = value of macro variable for replicate=0 ;
	tran_paramindex=4;
	varname="ssbs";
	run;
	
	/********************************************************************/
	/*** Dietary constituents #5 - vf: min. amount and lambda         ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var vf ;
	where vf > 0;
	output out=_min5(keep=minamount) min=minamount ;
	run;
	
	
	data _min5;
	set _min5;
	minamount=minamount/2;
	tran_paramindex=5;
	length varname $ 32;
	varname="vf";
	run;
	
	data _temp5;
	set inboxcox(keep=replicaterowid vf &covars.);
	if (vf > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp5, 
	 subject = replicaterowid, 
	 var     = vf,  
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
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=5;
	varname="vf";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #6 - pfab: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var pfab ;
	where pfab > 0;
	output out=_min6(keep=minamount) min=minamount ;
	run;
	
	
	data _min6;
	set _min6;
	minamount=minamount/2;
	tran_paramindex=6;
	length varname $ 32;
	varname="pfab";
	run;
	
	data _temp6;
	set inboxcox(keep=replicaterowid pfab &covars.);
	if (pfab > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp6, 
	 subject = replicaterowid, 
	 var     = pfab,  
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
	*tran_lambda=0.25;* = value of macro variable for replicate=0 ;
	tran_paramindex=6;
	varname="pfab";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #7 - other: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var other ;
	where other > 0;
	output out=_min7(keep=minamount) min=minamount ;
	run;
	
	
	data _min7;
	set _min7;
	minamount=minamount/2;
	tran_paramindex=7;
	length varname $ 32;
	varname="other";
	run;
	
	data _temp7;
	set inboxcox(keep=replicaterowid other &covars.);
	if (other > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp7, 
	 subject = replicaterowid, 
	 var     = other,  
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
	*tran_lambda=0.23;* = value of macro variable for replicate=0 ;
	tran_paramindex=7;
	varname="other";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #8 - water: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var water ;
	where water > 0;
	output out=_min8(keep=minamount) min=minamount ;
	run;
	
	
	data _min8;
	set _min8;
	minamount=minamount/2;
	tran_paramindex=8;
	length varname $ 32;
	varname="water";
	run;
	
	data _temp8;
	set inboxcox(keep=replicaterowid water &covars.);
	
	if (water > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp8, 
	 subject = replicaterowid, 
	 var     = water,  
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
	*tran_lambda=0.75;* = value of macro variable for replicate=0 ;
	tran_paramindex=8;
	varname="water";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #9 - milk: min. amount and lambda       ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var milk ;
	where milk > 0;
	output out=_min9(keep=minamount) min=minamount ;
	run;
	
	
	data _min9;
	set _min9;
	minamount=minamount/2;
	tran_paramindex=9;
	length varname $ 32;
	varname="milk";
	run;
	
	data _temp9;
	set inboxcox(keep=replicaterowid milk &covars.);
	if (milk > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp9, 
	 subject = replicaterowid, 
	 var     = milk,  
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
	*tran_lambda=0.34;* = value of macro variable for replicate=0 ;
	tran_paramindex=9;
	varname="milk";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #10 - mufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var mufa ;
	where mufa > 0;
	output out=_min10(keep=minamount) min=minamount ;
	run;
	
	data _min10;
	set _min10;
	minamount=minamount/2;
	tran_paramindex=10;
	length varname $ 32;
	varname="mufa";
	run;
	
	data _temp10;
	set inboxcox(keep=replicaterowid mufa &covars.);
	if (mufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp10, 
	 subject = replicaterowid, 
	 var     = mufa,  
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
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=10;
	varname="mufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #11 - pufa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var pufa ;
	where pufa > 0;
	output out=_min11(keep=minamount) min=minamount ;
	run;
	
	data _min11;
	set _min11;
	minamount=minamount/2;
	tran_paramindex=11;
	length varname $ 32;
	varname="pufa";
	run;
	
	data _temp11;
	set inboxcox(keep=replicaterowid pufa &covars.);
	if (pufa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp11, 
	 subject = replicaterowid, 
	 var     = pufa,  
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
	*tran_lambda=0.22;* = value of macro variable for replicate=0 ;
	tran_paramindex=11;
	varname="pufa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #12 - sfa: min. amount and lambda      ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var sfa ;
	where sfa > 0;
	output out=_min12(keep=minamount) min=minamount ;
	run;
	
	data _min12;
	set _min12;
	minamount=minamount/2;
	tran_paramindex=12;
	length varname $ 32;
	varname="sfa";
	run;
	
	data _temp12;
	set inboxcox(keep=replicaterowid sfa &covars.);
	
	if (sfa > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp12, 
	 subject = replicaterowid, 
	 var     = sfa,  
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
	*tran_lambda=0.33;* = value of macro variable for replicate=0 ;
	tran_paramindex=12;
	varname="sfa";
	run;
	
	
	/********************************************************************/
	/*** Dietary constituents #13 - sugars_free: min. amount and lambda ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var sugars_free ;
	where sugars_free > 0;
	output out=_min13(keep=minamount) min=minamount ;
	run;
	
	
	data _min13;
	set _min13;
	minamount=minamount/2;
	tran_paramindex=13;
	length varname $ 32;
	varname="sugars_free";
	run;
	
	data _temp13;
	set inboxcox(keep=replicaterowid sugars_free &covars.);
	if (sugars_free > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp13, 
	 subject = replicaterowid, 
	 var     = sugars_free,  
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
	*tran_lambda=0.44;* = value of macro variable for replicate=0 ;
	tran_paramindex=13;
	varname="sugars_free";
	run;


	/********************************************************************/
	/*** Dietary constituents #14 - energy: min. amount and lambda    ***/
	/********************************************************************/
	
	/* 1) Get minimum non-zero consumption amount */
	proc means data=datamrec min noprint;
	var energy ;
	where energy > 0;
	output out=_min14(keep=minamount) min=minamount ;
	run;
	
	data _min14;
	set _min14;
	minamount=minamount/2;
	tran_paramindex=14;
	length varname $ 32;
	varname="energy";
	run;
	
	data _temp14;
	set inboxcox(keep=replicaterowid energy &covars.);
	if (energy > 0) then output;
	run;
	
	/* 2) Call the <boxcox_survey> macro to find best normal transformation */
	options nonotes nosource;
	%boxcox_survey(
	 data    = _temp14, 
	 subject = replicaterowid, 
	 var     = energy,  
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
	*tran_lambda=0.22;* = value of macro variable for replicate=0 ;
	tran_paramindex=14;
	varname="energy";
	run;

 /*****************************************************************/
 /*** Final step: Append data sets of all dietary constituents  ***/
 /*****************************************************************/
	
	data work.xlambdas_rep&replicnum.;
	set _x1-_x14;
	run;
	
	data work.xminamount_rep&replicnum.;
	set _min1-_min14;
	label minamount= "Min. non zero amount divided by 2";
	run;


 /*************************************************************************/
 /*                                                                       */
 /* PRE NCI STEP: STANDARDIZE DATA <std_cov_boxcox24hr_conday_minamt>     */
 /*                                                                       */
 /*************************************************************************/

	/* Call the <std_cov_boxcox24hr_conday_minamt> macro with <xlambdas_rep&replicnum.> and <xminamount_rep&replicnum.> as inputs */
	%std_cov_boxcox24hr_conday_minamt(
	 data                       = datamrec, 
	 prestand_continuous_covars = ageat24hr deprivation sedentary0 bmi0 urine_na_24h , 
	 rec24hr_epis_vars          = pfpb rg wg ssbs, 
	 rec24hr_daily_vars         = vf pfab other water milk mufa pufa sfa sugars_free energy, 
	 boxcox_tran_lambda_data    = xlambdas_rep&replicnum., 
	 minamount_data             = xminamount_rep&replicnum., 
	 print                      = y, 
	 titles                     = 3 ); 
	 

	/* note: contivuous covariates in <prestand_continuous_covars> are standardized and renamed as std_<covar> */

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 1: Fit the multivariate meas. error model (MULTIVAR_MCMC)    */
 /*                                                                       */
 /*************************************************************************/

	/* note: multivar_mcmc may take 24 to 48h to complete */

	%multivar_mcmc(
	 data                        = stdcov_stdbc24hr_conday_out, 
	 subject                     = replicaterowid, 
	 weight_var                  =  , 
	 repeat                      = r24_no, 
	 conday_epis_vars            = conday_pfpb  conday_rg  conday_wg  conday_ssbs, 
	 gst_rec24hr_epis_vars       = stdbc_pfpb  stdbc_rg  stdbc_wg  stdbc_ssbs, 
	 gst_rec24hr_daily_vars      = stdbc_vf  stdbc_pfab  stdbc_other  stdbc_water  stdbc_milk  stdbc_mufa  stdbc_pufa  stdbc_sfa  
	stdbc_sugars_free  stdbc_energy, 
	 covars_epis_prob            = constant1   sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  famhist0_1 
	famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  riskfactor_1 
	riskfactor_2  physact0_2 physact0_3 std_ageat24hr  std_deprivation  std_sedentary0  std_bmi0  std_urine_na_24h , 
	 covars_epis_amt             = constant1   sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  famhist0_1 
	famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  riskfactor_1 
	riskfactor_2  physact0_2 physact0_3 std_ageat24hr  std_deprivation  std_sedentary0  std_bmi0  std_urine_na_24h , 
	 covars_daily_amt            = constant1   sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  famhist0_1 
	famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  riskfactor_1 
	riskfactor_2  physact0_2 physact0_3 std_ageat24hr  std_deprivation  std_sedentary0  std_bmi0  std_urine_na_24h , 
	 nev_consumers_epis1         = , 
	 covars_prob_consumer_epis1  = , 
	 set_seed_mcmc               = 85600495, 
	 set_number_mcmc_iterations  = 2500, 
	 set_number_burn_iterations  = 500, 
	 set_thin                    = 10, 
	 prior_sigmau_mean_data      = , 
	 sigmau_constant             = , 
	 gen_inverse                 = , 
	 print                       = y, 
	 titles                      = 3, 
	 std_print_store             = y, 
	 notes_print                 = y, 
	 out_lib                     = &outlib., 
	 out_store_label             = mcmc_rep&replicnum., 
	 out_save_label_max5char     = _r&replicnum., 
	 set_number_saved_out_data   = , 
	 save_mcmc_u_out_data        = , 
	 set_number_post_mcmc_u_out  = 1000,
	 traceplots_method1_gpath    = , 
	 traceplots_method2_file_pdf = epi_trace_rep&replicnum..pdf, 
	 optional_iml_store_data     = backtran_out, 
	 optional_iml_store_names    = constant1   sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  famhist0_1 
	famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  riskfactor_1 
	riskfactor_2  physact0_2 physact0_3 std_ageat24hr  std_deprivation  std_sedentary0  std_bmi0  std_urine_na_24h 
	tran_paramindex tran_lambda tran_center tran_scale minamount 
	 ); 
	
	
	/* Save transformations and min. amount of dietary constituents */
	data &outlib..backtran_out&replicnum.;
	set work.backtran_out;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* NCI PART 2: Simulation of pseudo-individuals (MULTIVAR_DISTRIB)       */
 /*                                                                       */
 /*************************************************************************/

	/* note: the simulation data is divided in 5 subsets to lower
		memory requirements (200 simulations per subset). The number of simulations must
		be consistent with the value of <set_number_post_mcmc_u_out> in
		the multivar_mcmc macro (i.e., 1000 in the present code). */
	
	%let num_subsets_monte_carlo_distrib = 5;
    %let subset_number_post_mcmc_u_out = 200;

	/* begin subset loop */
	%do subset_monte_carlo_distrib = 1 %to &num_subsets_monte_carlo_distrib;

	/* create seed for current subset */
	%let loopseed_distrib = %eval(89009890 + (&replicnum * 10000) + (&subset_monte_carlo_distrib * 111));
	
	/*  Handling subset # <subset_monte_carlo_distrib> / 5 */
	data &outlib..subset_r&replicnum;
		set &outlib..multivar_post_mcmc_u_out_r&replicnum;
	* keep only pseudo-individuals of current subset;
	if (&subset_number_post_mcmc_u_out * (&subset_monte_carlo_distrib - 1)) < iteration <= (&subset_number_post_mcmc_u_out * &subset_monte_carlo_distrib);
	iteration = iteration - (&subset_number_post_mcmc_u_out * (&subset_monte_carlo_distrib - 1));
	run;
	 
	 /* Call the <multivar_distrib> macro to predict usual intakes for current subset*/
	 %multivar_distrib(
	 multivar_mcmc_out_lib           = &outlib. ,  
	 multivar_mcmc_out_store_label   = mcmc_rep&replicnum., 
	 t_weightavg_covariates_list1    = constant1   sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  famhist0_1 
	famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  riskfactor_1 
	riskfactor_2  physact0_2 physact0_3 std_ageat24hr  std_deprivation  std_sedentary0  std_bmi0  std_urine_na_24h ,  
	 t_weightavg_covariates_list2    = , 
	 set_value_for_weight_cov_list1  = , 
	 set_value_for_weight_cov_list2  = , 
	 nev_consumers_epis1             = , 
	 covars_prob_consumer_epis1      = , 
	 optional_input_data             = subj1rec , 
	 optional_input_data_var_list    = replicaterowid constant1 sexmeno_2 sexmeno_3 sexmeno_4  region2_Scotland region2_Wales  
	famhist0_1 famhist0_2   employ0_2 employ0_3  university_1  smk0_1 smk0_2  roh0_1  rx0_1 rx0_2  hrt0_1  dietsuppl_1  
	riskfactor_1 riskfactor_2  physact0_2 physact0_3 ageat24hr deprivation sedentary0 bmi0 urine_na_24h urine_na_24h  , 
	 optional_input_mcmc_u_out_data  = subset_r&replicnum., 
	 additional_output_var_list      = urine_na_24h , 
	 additional_output_subject_var   = replicaterowid , 
	 output_mcmc_weight_var          = y  , 
	 set_seed_distrib                = &loopseed_distrib, 
	 set_number_monte_carlo_rand_obs = ,  
	 print                           = y 
	 ); 
	
	 /*Save current subset (i.e., pseudo-individuals data) */
	data &outlib..allsubsets_mc_t_distrib&replicnum.;
	%if (&subset_monte_carlo_distrib = 1) %then %do;
		set mc_t_distrib_out;
	%end;
 	%else %do;
		set &outlib..allsubsets_mc_t_distrib&replicnum. mc_t_distrib_out;
	%end;
	run;

    %end; /* end of monte carlo simulation subset loop */

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: Format output data from MULTIVAR_DISTRIB               */
 /*                                                                       */
 /*************************************************************************/

	data &outlib..usintake_mc_t&replicnum.(drop=mc_t1-mc_t14 i mc_prob: mc_amount: );
	set &outlib..allsubsets_mc_t_distrib&replicnum.;
	* Assign original variable names to the output from the MULTIVAR_DISTRIB macro ;
	array outmc (*) mc_t1-mc_t14 ;
	array clean (*) pfpb rg wg ssbs vf pfab other water milk mufa pufa sfa sugars_free energy ;
	do i=1 to dim(outmc);
		clean(i)=outmc(i);
	end;
	run;
	
	proc sort data=&outlib..usintake_mc_t&replicnum.;
		by replicaterowid;
	run;

	/* clean temporary data */
	proc datasets lib=work nolist nodetails;
	delete mc_t_distrib_out ;
	run;

	proc datasets lib=&outlib. nolist nodetails;
	delete multivar_post_mcmc_u_out_r&replicnum subset_r&replicnum allsubsets_mc_t_distrib&replicnum.;
	run;
	
 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: Apply the HEFI-2019 scoring algorithm                  */
 /*                                                                       */
 /*************************************************************************/

 /* Call the <HEFI2019> macro to apply the HEFI-2019 scoring algorithm*/
	%HEFI2019(
	 indata             = &outlib..usintake_mc_t&replicnum., 
	 vegfruits          = vf, 
	 wholegrfoods       = wg, 
	 nonwholegrfoods    = rg, 
	 profoodsanimal     = pfab, 
	 profoodsplant      = pfpb, 
	 otherfoods         = other, 
	 waterhealthybev    = water, 
	 unsweetmilk        = milk, 
	 unsweetplantbevpro = 0, 
	 otherbeverages     = ssbs, 
	 mufat              = mufa, 
	 pufat              = pufa, 
	 satfat             = sfa, 
	 freesugars         = sugars_free, 
	 sodium             = urine_na_24h, 
	 energy             = energy, 
	 outdata            = &outlib..usintake_mc_t&replicnum. 
	 ); 


 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: Restricted cubic spline transformation                 */
 /*                                                                       */
 /*************************************************************************/

	/* From distribution analysis based on usual intakes
	   # HEFI2019_TOTAL_SCORE .05 spline knot 1 = 29.27 
	   # HEFI2019_TOTAL_SCORE .35 spline knot 2 = 42.53
	   # HEFI2019_TOTAL_SCORE .65 spline knot 3 = 50.35
	   # HEFI2019_TOTAL_SCORE .95 spline knot 4 = 60.70 */

	/* # ENERGY .05 spline knot 1 = 1.46
	   # ENERGY .35 spline knot 2 = 1.92 
	   # ENERGY .65 spline knot 3 = 2.26
	   # ENERGY .95 spline knot 4 = 2.95 */


	data &outlib..usintake_mc_t&replicnum. ;
	set &outlib..usintake_mc_t&replicnum.  ;
	
	************************;
	* Total HEFI-2019 score ;
	************************;
	hefi2019_rcs_lin = HEFI2019_TOTAL_SCORE;
	
	hefi2019_rcs_s1  = (HEFI2019_TOTAL_SCORE > 29.27)*(HEFI2019_TOTAL_SCORE - 29.27)**3 - (HEFI2019_TOTAL_SCORE > 50.35) * (( 60.7 - 
	29.27)/( 60.7 - 50.35)) * (HEFI2019_TOTAL_SCORE - 50.35)**3 + (HEFI2019_TOTAL_SCORE > 60.7) * (( 50.35 - 29.27)/( 60.7 - 50.35)) * 
	(HEFI2019_TOTAL_SCORE - 60.7)**3 ;
	
	hefi2019_rcs_s2  = (HEFI2019_TOTAL_SCORE > 42.53)*(HEFI2019_TOTAL_SCORE - 42.53)**3 - (HEFI2019_TOTAL_SCORE > 50.35) * (( 60.7 - 
	42.53)/( 60.7 - 50.35)) * (HEFI2019_TOTAL_SCORE - 50.35)**3 + (HEFI2019_TOTAL_SCORE > 60.7) * (( 50.35 - 42.53)/( 60.7 - 50.35)) * 
	(HEFI2019_TOTAL_SCORE - 60.7)**3 ;
	
	************************;
	* Total energy intakes ;
	************************;
	* rescale energy to 1000 calories;
	ENERGY = ENERGY/1000;
	
	energy_rcs_lin = ENERGY;
	
	energy_rcs_s1 = (ENERGY > 1.46)*(ENERGY - 1.46)**3 - (ENERGY > 2.26) * (( 2.95 - 1.46)/( 2.95 - 2.26)) * 
	(ENERGY - 2.26)**3 + (ENERGY > 2.95) * (( 2.26 - 1.46)/( 2.95 - 2.26)) * (ENERGY - 2.95)**3 ;
	
	energy_rcs_s2 = (ENERGY > 1.92)*(ENERGY - 1.92)**3 - (ENERGY > 2.26) * (( 2.95 - 1.92)/( 2.95 - 2.26)) * 
	(ENERGY - 2.26)**3 + (ENERGY > 2.95) * (( 2.26 - 1.92)/( 2.95 - 2.26)) * (ENERGY - 2.95)**3 ;
	
	run;

	/* Restricted cubic spline transformation code from: 
	Desquilbet, Loic, and Francois Mariotti. Dose-Response Analyses Using Restricted Cubic
	Spline Functions in Public Health Research. Statistics in Medicine, 2010.
	https://doi.org/10.1002/sim.3841.*/

 /*************************************************************************/
 /*                                                                       */
 /* POST NCI STEP: Average intakes across pseudo-individuals              */
 /*                                                                       */
 /*************************************************************************/

	proc means data=&outlib..usintake_mc_t&replicnum. noprint;
	by replicaterowid;
	var pfpb rg wg ssbs vf pfab other water milk mufa pufa sfa sugars_free energy
	RATIO_VF RATIO_WGTOT RATIO_WGGR RATIO_PRO RATIO_PLANT RATIO_FA RATIO_BEV SFA_PERC SUG_PERC SODDEN
	HEFI2019C1_VF HEFI2019C2_WHOLEGR HEFI2019C3_GRRATIO HEFI2019C4_PROFOODS HEFI2019C5_PLANTPRO
	HEFI2019C6_BEVERAGES HEFI2019C7_FATTYACID HEFI2019C8_SFAT HEFI2019C9_FREESUGARS HEFI2019C10_SODIUM
	HEFI2019_TOTAL_SCORE 
	HEFI2019_RCS_lin HEFI2019_RCS_S1 HEFI2019_RCS_S2 
	ENERGY_RCS_lin ENERGY_RCS_S1 ENERGY_RCS_S2 ;
	output out=&outlib..usintake_mc_t&replicnum. mean= ;
	run;

	/* merge predicted usual intakes and derived variables with 1-row per subject data */
	data &outlib..subj1recres&replicnum.;
	merge subj1rec
		  &outlib..usintake_mc_t&replicnum.(keep = replicaterowid pfpb rg wg ssbs vf pfab other water milk mufa pufa sfa sugars_free energy
			RATIO_VF RATIO_WGTOT RATIO_WGGR RATIO_PRO RATIO_PLANT RATIO_FA RATIO_BEV SFA_PERC SUG_PERC SODDEN
			HEFI2019C1_VF HEFI2019C2_WHOLEGR HEFI2019C3_GRRATIO HEFI2019C4_PROFOODS HEFI2019C5_PLANTPRO
			HEFI2019C6_BEVERAGES HEFI2019C7_FATTYACID HEFI2019C8_SFAT HEFI2019C9_FREESUGARS HEFI2019C10_SODIUM
			HEFI2019_TOTAL_SCORE 
			HEFI2019_RCS_lin HEFI2019_RCS_S1 HEFI2019_RCS_S2 
			ENERGY_RCS_lin ENERGY_RCS_S1 ENERGY_RCS_S2 );
	by replicaterowid;
	run;
	
	/* Close log */
	proc printto; run;
	
	%end; /* end of bootstrap loop */
	
	%mend replicateloop;

	%replicateloop(replicfirst = 0,
				   repliclast  = 0 );
				   
	/* note: for simplicity, the present example is for the original sample (replicate = 0) */


	 /*************************************************************************/
	 /*                                                                       */
 	 /*                          BOOTSTRAP LOOP END                           */
	 /*                                                                       */
	 /*************************************************************************/


 /*************************************************************************/
 /*                                                                       */
 /* Append original and bootstrap `subj1recres` data sets                 */
 /*                                                                       */
 /*************************************************************************/

	/* note: code below assumes <replicateloop> has been executed for all 250 bootstrap samples */

/* Append base and bootstrap subj1recres data */
	data temp.subj1recres ;
		set baselib.subj1recres0
			bootlib.subj1recres1-bootlib.subj1recres250
			;
			
	* rescale energy back to calories ;
	energy=energy*1000;

	drop
	/* calibration modeling variables */
		constant1
	/* intakes on a given day - not needed*/
		vf1-vf5 pfab1-pfab5 other1-other5 water1-water5 milk1-milk5
		mufa1-mufa5 pufa1-pufa5 sfa1-sfa5 sugars_free1-sugars_free5 energy1-energy5 
		pfpb1-pfpb5 rg1-rg5  wg1-wg5 ssbs1-ssbs5 ;
	run;

	proc sort data=temp.subj1recres ;
		by subjectid ; /* original participant identifier */
	run;
	
	/* Add outcome data to mesurement error-corrected data */
	data fmtdata.postNCI;
		retain replicate replicaterowid subjectid;
		merge fmtdata.screening(keep=subjectid time_to_primary event_primary 
				cens_death_primary rename=(nomat=subjectid))
			  fmtdata.lost_follow_up(keep=nomat lost_flag  rename=(nomat=subjectid))
			  temp.subj1recres(in=nci);
		by subjectid;
	* change <lost_flag> value for participants that were NOT lost to follow-up ;
		if missing(lost_flag) then lost_flag=0;
	* recode deaths AND lost to follow-up under one common variable;
		cens_primary = cens_death_primary;
		if lost_flag=1 then cens_primary = 1;
	* drop old censoring variables ;
		drop cens_death_primary lost_flag;
	* keep only included participants;
	if nci then output;
	run;
	
	proc sort data=fmtdata.postNCI;
		by replicate replicaterowid;
	run;
	

 /*************************************************************************/
 /*                                                                       */
 /*                              END OF CODE                              */
 /*                                                                       */
 /*************************************************************************/