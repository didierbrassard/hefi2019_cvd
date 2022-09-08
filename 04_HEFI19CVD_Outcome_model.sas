 /*************************************************************************/
 /*                                                                       */
 /*               HEFI-2019 and CVD risk in the UK Biobank                */
 /*                                                                       */
 /*                  Brassard et al. Am J Clin Nutr 2022                  */
 /*               Code 4: Outcome model (Cox PH regression)               */
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
	%let suffix=_f ;
	
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
 /* HEFI19: OUTPUT DISTRIBUTION PERCENTILES AND RCS KNOTS (NCI)           */
 /*                                                                       */
 /*************************************************************************/

/* HEFI-2019: output knots value used to RCS-transform data (knots were estimated in the overall sample) */
	proc transpose data=nci.hefi_spl(where=(food="HEFI2019_TOTAL_SCORE")) 
			out=knots_rcs_hefi(rename=(col1=value));
		var k1-k4;
	run;
	
	/* HEFI-2019: Show knots value (i.e., percentile of the distribution used to derive the restricted cubic spline) */
	proc print data=knots_rcs_hefi;run;
	
/* HEFI-2019: transpose predetermined percentile values for survival curves */
	proc transpose data=distrib.distribtotal_w_cens_death0(Where=(cens_death=0)) 
			out=_contrast_cens_death0_t(rename=(col1=hefi2019_total_score)) name=pctile;
		var Pctile5 Pctile10 Pctile25 Pctile50 Pctile75 Pctile90 Pctile95;
	run;
	
	/* HEFI-2019: show percentile values */
	proc print data=_contrast_cens_death0_t;run;
	
	/************************************************/
	/*  Transform `raw` hefi2019 to `rcs`-modelled  */
	/************************************************/
	
	/* note: this step is needed since exposure is restricted cubic spline (RCS)-transformed */
			
	data x_to_splx;
		hefi2019 = 29.3 ; * 5th percentile ;
		output;
		hefi2019 = 33.1 ; * 10th percentile ;
		output;
		hefi2019 = 39.5 ; * 25th percentile ;
		output;
		hefi2019 = 46.6 ; * 50th percentile ;
		output;
		hefi2019 = 53.1 ; * 75th percentile ;
		output;
		hefi2019 = 58.1 ; * 90th percentile ;
		output;
		hefi2019 = 60.7 ; * 95th percentile ;
		output;
	run;

	/* From distribution analysis based on usual intakes: */
	/* # HEFI2019_TOTAL_SCORE .05 spline knot 1 = 29.27 */
	/* # HEFI2019_TOTAL_SCORE .35 spline knot 2 = 42.53 */
	/* # HEFI2019_TOTAL_SCORE .65 spline knot 3 = 50.35 */
	/* # HEFI2019_TOTAL_SCORE .95 spline knot 4 = 60.70 */

	data rcs_hefi;
	set x_to_splx;
	
	hefi2019_RCS_lin = hefi2019;
	
	hefi2019_RCS_S1 = (hefi2019 > 29.27)*(hefi2019 - 29.27)**3 - (hefi2019 > 50.35) * ((60.70 - 29.27)/(60.70 - 
	 50.35)) * (hefi2019 - 50.35)**3 + (hefi2019 > 60.70) * ((50.35 - 29.27)/(60.70 - 50.35)) * (hefi2019 - 60.70)**3 ;
	
	hefi2019_RCS_S2 = (hefi2019 > 42.53)*(hefi2019 - 42.53)**3 - (hefi2019 > 50.35) * ((60.70 - 42.53)/(60.70 - 
	 50.35)) * (hefi2019 - 50.35)**3 + (hefi2019 > 60.70) * ((50.35 - 42.53)/(60.70 - 50.35)) * (hefi2019 - 60.70)**3 ;
	run;

	/* Restricted cubic spline transformation code from: 
	Desquilbet, Loic, and Francois Mariotti. Dose-Response Analyses Using Restricted Cubic
	Spline Functions in Public Health Research. Statistics in Medicine, 2010.
	https://doi.org/10.1002/sim.3841.*/


 /*************************************************************************/
 /*                                                                       */
 /* ENERGY: OUTPUT MEAN AND RCS KNOTS (BASED ON NCI DISTRIBUTION)         */
 /*                                                                       */
 /*************************************************************************/

/* ENERGY: output knots value used to RCS-transform data (knots were estimated in the overall sample) */
	proc transpose data=nci.hefi_spl(where=(food="energy")) 
			out=knots_rcs_energy(rename=(col1=value));
		var k1-k4;
	run;

	/* ENERGY: Show knots value (i.e., percentile of the distribution used to derive the restricted cubic spline) */
	proc print data=knots_rcs_energy;run;

/* ENERGY: output mean energy */
	data _energy;
		set distrib.distribraw_w0(where=(varname="energy"));
	* values for energy were divided /1000;
	Mean = round(Mean,100)/1000 ;
	run;

	/* ENERGY: show mean */
	proc print data=_energy;run;
	
/* ENERGY: convert `raw` mean energy to `rcs`-modeled values */
	
	/* note: for energy, we use the same mean intake for all predetermined percentile of hefi2019 score */

	data e_to_sple;
	energy = 2.1 ;
	run;

	/* From distribution analysis based on usual intakes: */
	/* # ENERGY .05 spline knot 1 = 1.46 */
	/* # ENERGY .35 spline knot 2 = 1.92 */
	/* # ENERGY .65 spline knot 3 = 2.26 */
	/* # ENERGY .95 spline knot 4 = 2.95 */

	data rcs_e;
	set e_to_sple;
	
	energy_RCS_lin = energy;
	
	energy_RCS_S1 = (energy > 1.46)*(energy - 1.46)**3 - (energy > 2.26) * ((2.95 - 1.46)/(2.95 - 2.26)) * 
	 (energy - 2.26)**3 + (energy > 2.95) * ((2.26 - 1.46)/(2.95 - 2.26)) * (energy - 2.95)**3 ;
	
	energy_RCS_S2 = (energy > 1.92)*(energy - 1.92)**3 - (energy > 2.26) * ((2.95 - 1.92)/(2.95 - 2.26)) * 
	 (energy - 2.26)**3 + (energy > 2.95) * ((2.26 - 1.92)/(2.95 - 2.26)) * (energy - 2.95)**3 ;
	run;
	
	/* Restricted cubic spline transformation code from: 
	Desquilbet, Loic, and Francois Mariotti. Dose-Response Analyses Using Restricted Cubic
	Spline Functions in Public Health Research. Statistics in Medicine, 2010.
	https://doi.org/10.1002/sim.3841.*/
	
 /*************************************************************************/
 /*                                                                       */
 /* COMBINE HEFI2019 AND ENERGY EXPOSURE DATA FOR OUTCOME MODEL           */
 /*                                                                       */
 /*************************************************************************/
	
	data InExposure reslib.InExposure ;
	retain p pctile;
	if _N_=1 then set rcs_e; * keep energy constant across all predetermined percentiles;
		merge rcs_hefi _contrast_cens_death0_t(keep=pctile);
	* identify scores based on their percentile;
	p=input(compress(pctile,,'a'),10.);
	run;
	
	/* note: InExposure data is used to output survival curves
		at predetermined percentiles of hefi2019 score */

 /*************************************************************************/
 /*                                                                       */
 /* OUTCOME MODEL: WEIGHTED (SW) COX PH REG WITH RCS HEFI19 AND ENERGY    */
 /*                                                                       */
 /*************************************************************************/

	proc phreg data=fmtdata.sw_final outest=reslib.coxparm_w;
		by replicate; * for bootstrap variance estimation ;
		model time_to_primary * event_primary(0) = hefi2019_RCS_lin hefi2019_RCS_S1 hefi2019_RCS_S2
			energy_RCS_lin energy_RCS_S1 energy_RCS_S2 / ties=efron;
		weight sw;
		baseline covariates=InExposure survival=Survival out=SurvivalPlot / rowid=pctile;
		* note: <baseline> statement used to output survival curves at predetermined percentiles
			and constant energy intake, based on Cox model parameters estimated in 
			sw-weighted sample (pseudopopulation) ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* SURVIVAL CURVES: SAVE PR(SURVIVAL) AND PERFORM SOME DATA MANIPULATION */
 /*                                                                       */
 /*************************************************************************/

/* save survival curves (tall/long data, without s(t) contrast) */
	data reslib.survivalplot_t;
		set SurvivalPlot;
	* make numerical percentile;
	p=input(compress(pctile,,'a'),10.);
	run;
	
	proc sort data=reslib.survivalplot_t;
		by replicate p time_to_primary;
	run;

/* save survival curve (wide data, without s(t) contrast) */
	proc transpose data=reslib.survivalplot_t(where=(survival ne 1)) 
			out=reslib.survivalplot_w1(drop=_NAME_) prefix=t;
		by replicate p;
		var survival;
		id time_to_primary;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* SURVIVAL CURVES: OUTPUT EXPOSURE CONTRASTS BASED ON PERCENTILES       */
 /*                                                                       */
 /*************************************************************************/
		
/* transpose s(t) for each exposure level (i.e., hefi percentile) */
	proc sort data=reslib.survivalplot_t out=_insurvival;
	by replicate time_to_primary ;
	run;
	 
	proc transpose data=_insurvival out=st_exp1(rename=(p1=p5)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=5);
	run;
	 
	proc transpose data=_insurvival out=st_exp2(rename=(p1=p10)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=10);
	run;
	 
	proc transpose data=_insurvival out=st_exp3(rename=(p1=p25)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=25);
	run;
	 
	proc transpose data=_insurvival out=st_exp4(rename=(p1=p50)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=50);
	run;
	 
	proc transpose data=_insurvival out=st_exp5(rename=(p1=p75)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=75);
	run;
	 
	proc transpose data=_insurvival out=st_exp6(rename=(p1=p90)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=90);
	run;
	 
	proc transpose data=_insurvival out=st_exp7(rename=(p1=p95)) prefix=p;
	by replicate time_to_primary ;
	var survival ;
	where (survival ne 1) & (p=95);
	run;

/* merge all (wide) exposure data, calculate risk difference and ratio */
	data reslib.survivalplot_w2(drop=i rename=(_NAME_=name _LABEL_=label));
	retain replicate ;
	merge st_exp1-st_exp7;
	by replicate time_to_primary;
	
	* change Pr(survival) to Pr(outcome); 
	array p(*) p5 p10 p25 p50 p75 p90 p95 ;
	do i=1 to dim(p);
		p(i) = 1-p(i);
	end;
	
	* Operation 1: risk ratio for RR_p50_p5 ;
	 RR_p50_p5 = (p50/p5);
	* Operation 1: risk difference for RD_p50_p5 ;
	 RD_p50_p5 = (p50-p5);
	 
	* Operation 2: risk ratio for RR_p50_p10 ;
	 RR_p50_p10 = (p50/p10);
	* Operation 2: risk difference for RD_p50_p10 ;
	 RD_p50_p10 = (p50-p10);
	 
	* Operation 3: risk ratio for RR_p50_p25 ;
	 RR_p50_p25 = (p50/p25);
	* Operation 3: risk difference for RD_p50_p25 ;
	 RD_p50_p25 = (p50-p25);
	 
	* Operation 4: risk ratio for RR_p75_p50 ;
	 RR_p75_p50 = (p75/p50);
	* Operation 4: risk difference for RD_p75_p50 ;
	 RD_p75_p50 = (p75-p50);
	 
	* Operation 5: risk ratio for RR_p90_p50 ;
	 RR_p90_p50 = (p90/p50);
	* Operation 5: risk difference for RD_p90_p50 ;
	 RD_p90_p50 = (p90-p50);
	 
	* Operation 6: risk ratio for RR_p95_p50 ;
	 RR_p95_p50 = (p95/p50);
	* Operation 6: risk difference for RD_p95_p50 ;
	 RD_p95_p50 = (p95-p50);
	 
	label _NAME_="Outcome";
	run;

 /*************************************************************************/
 /*                                                                       */
 /*                              END OF CODE                              */
 /*                                                                       */
 /*************************************************************************/