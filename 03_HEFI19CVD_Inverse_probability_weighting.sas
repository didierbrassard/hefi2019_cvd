 /*************************************************************************/
 /*                                                                       */
 /*               HEFI-2019 and CVD risk in the UK Biobank                */
 /*                                                                       */
 /*                  Brassard et al. Am J Clin Nutr 2022                  */
 /*                 Code 3: Inverse probability weighting                 */
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
 /* IPW_CONT_NORMAL: DENOMINATOR MODEL (X | L )                           */
 /*                                                                       */
 /*************************************************************************/

	/* note: proc surveyreg is used to transform continuous variables to RCS
		with the <effect> statement, but other procedures could also be used for
		standard-normal IPW estimation */

	ods select none;
	
	proc surveyreg data=fmtdata.postNCI;
		by replicate;
		effect spl_1=spline(ageat24hr / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 27.5 50 72.5 95));
		effect spl_2=spline(deprivation / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(10 50 90));
		effect spl_3=spline(sedentary0 / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 35 65 95));
		effect spl_4=spline(bmi0 / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 27.5 50 72.5 95));
		model hefi2019_total_score = sexmeno_2 sexmeno_3 sexmeno_4 region2_Scotland 
			region2_Wales famhist0_1 famhist0_2 employ0_2 employ0_3 university_1 smk0_1 
			smk0_2 roh0_1 rx0_1 rx0_2 hrt0_1 dietsuppl_1 riskfactor_1 riskfactor_2 
			physact0_2 physact0_3 spl_1 spl_2 spl_3 spl_4 / adjrsq;
		ods output FitStatistics=_ss_d;
		output out=_den predicted=xb_d;
	run;
	
	ods select all;
	
 /*************************************************************************/
 /*                                                                       */
 /* IPW_CONT_NORMAL: DENOMINATOR PROBABILITY DENSITY FUNCTION             */
 /*                                                                       */
 /*************************************************************************/

	data _ss_d(keep=replicate nValue1 rename=(nValue1 = root_mse_d));
		set _ss_d (where=(Label1 in ("Root MSE" "Racine MSE")));
	run;

	data est_dens_d;
	merge _den _ss_d;
	by replicate ;
	pdf_den_denom = pdf("normal",hefi2019_total_score,xb_d,root_mse_d);
	drop xb_d root_mse_d ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* IPW_CONT_NORMAL: NUMERATOR MODEL (X | __ )                            */
 /*                                                                       */
 /*************************************************************************/

	ods select none; 
	 
	proc surveyreg data=fmtdata.postNCI; 
		by replicate; 
		model hefi2019_total_score=; 
		ods output FitStatistics=_ss_n(where=(Label1 in ("Root MSE" "Racine MSE"))); 
		output out=_num predicted=xb_n; 
	run; 
	 
	ods select all;
	
 /*************************************************************************/
 /*                                                                       */
 /* IPW_CONT_NORMAL: NUMERATOR PROBABILITY DENSITY FUNCTION               */
 /*                                                                       */
 /*************************************************************************/

	data est_dens_n(keep=replicate replicaterowid pdf_den_numer );
	merge _num _ss_n(keep=replicate nValue1 rename=(nValue1=root_mse_n));
	by replicate;
	pdf_den_numer = pdf("normal",hefi2019_total_score,xb_n,root_mse_n);
	run;
	
 /*************************************************************************/
 /*                                                                       */
 /* IPW_CONT_NORMAL: WEIGHT ESTIMATION                                    */
 /*                                                                       */
 /*************************************************************************/

	proc sort data=est_dens_d;
		by replicate replicaterowid;
	run;
	
	proc sort data=est_dens_n;
		by replicate replicaterowid;
	run;
	
	data fmtdata.sw_cont;
	merge est_dens_d est_dens_n ;
	by replicate replicaterowid;
	* calculate stabilized weight based on numerator and denominator ;
	sw_a= pdf_den_numer / pdf_den_denom;
	label sw_a = "Inverse-probability-weighting for hefi2019_total_score" ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* IPW_CENSORING: DENOMINATOR MODEL for censoring (C | X, E, L)          */
 /*                                                                       */
 /*************************************************************************/

	ods select none;
	
	proc logistic data=fmtdata.sw_cont;
		by replicate;
		effect spl_1=spline(ageat24hr / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 27.5 50 72.5 95));
		effect spl_2=spline(deprivation / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(10 50 90));
		effect spl_3=spline(sedentary0 / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 35 65 95));
		effect spl_4=spline(bmi0 / details naturalcubic basis=tpf(noint) 
			knotmethod=percentilelist(5 35 65 95));
		model cens_primary = hefi2019_rcs_lin hefi2019_rcs_s1 hefi2019_rcs_s2 
			energy_rcs_lin energy_rcs_s1 energy_rcs_s2 sexmeno_2 sexmeno_3 sexmeno_4 
			region2_Scotland region2_Wales famhist0_1 famhist0_2 employ0_2 employ0_3 
			university_1 smk0_1 smk0_2 roh0_1 rx0_1 rx0_2 hrt0_1 dietsuppl_1 riskfactor_1 
			riskfactor_2 physact0_2 physact0_3 spl_1 spl_2 spl_3 spl_4 / expb;
		output out=est_prob_d_cens (keep=replicate replicaterowid pd_cens) p=pd_cens;
	run;
	
	ods select all;
	
	proc sort data=est_prob_d_cens;
	by replicate replicaterowid ;
	run;


 /*************************************************************************/
 /*                                                                       */
 /* IPW_CENSORING: NUMERATOR MODEL (C | X)                                */
 /*                                                                       */
 /*************************************************************************/
	
	ods select none;
	
	proc logistic data=fmtdata.sw_cont;
		by replicate;
		ods exclude ClassLevelInfo ModelAnova Association GlobalTests Oddsratios;
		model cens_primary=hefi2019_rcs_lin hefi2019_rcs_s1 hefi2019_rcs_s2;
		output out=est_prob_n_cens p=pn_cens;
	run;
	
	ods select all;
	
	proc sort data=est_prob_n_cens;
		by replicate replicaterowid;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* IPW_CENSORING: WEIGHT ESTIMATION                                      */
 /*                                                                       */
 /*************************************************************************/

	data fmtdata.sw_cens;
		merge est_prob_d_cens est_prob_n_cens;
		by replicate replicaterowid ;
		* calculate stabilized IPCW weight ;
		sw_c= pn_cens / pd_cens;
		label sw_c = "Inverse-probability-of-censoring-weight. (adj. for subj. with cens_primary=1)";
		if (cens_primary=0) then output;
		* note: observations with censoring=1 only contribute to censoring models;
	run;

 /*************************************************************************/
 /*                                                                       */
 /* CALCULATE SW_A*SW_C PRODUCT TO OBTAIN COMBINED WEIGHT                 */
 /*                                                                       */
 /*************************************************************************/

	data fmtdata.sw_final;
		set fmtdata.sw_cens;
	* make final set of weights based on sw_a and sw_c;
		sw = sw_a * sw_c;
	label sw = "Total HEFI-2019 inverse probability weight (IPTW*IPCW)" ;
	run;

	proc sort data=fmtdata.sw_final;
		by replicate replicaterowid;
	run;
	

 /*************************************************************************/
 /*                                                                       */
 /*                              END OF CODE                              */
 /*                                                                       */
 /*************************************************************************/
