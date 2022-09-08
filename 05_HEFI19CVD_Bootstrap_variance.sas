 /*************************************************************************/
 /*                                                                       */
 /*               HEFI-2019 and CVD risk in the UK Biobank                */
 /*                                                                       */
 /*                  Brassard et al. Am J Clin Nutr 2022                  */
 /*                 Code 5: Bootstrap variance estimation                 */
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
 /* SURVIVAL CURVES: BOOTSTRAP VARIANCE ESTIMATION                        */
 /*                                                                       */
 /*************************************************************************/

/* separate original and bootstrap replicates */
	data survivalplot0 survivalplotbs;
		set reslib.survivalplot_w1;
	if replicate=0 then output survivalplot0;
		else output survivalplotbs;
	run;

	proc sort data=survivalplotbs;
		by p;
	run;
	 
	proc sort data=survivalplot0;
		by p;
	run;

/* output number of successful (non-missing) bootstraps */
	proc means data=survivalplotbs nmiss noprint;
	by p ;
	var t: ;
	output out=_nboot(drop=_TYPE_) nmiss=;
	run;
	
	/* of note, all time values are coded as t<months>,
		e.g., t12 refers to survival at 12 months */
	 
	data _nboot;
	set _nboot;
	array varlist (*) t: ;
	* Number of non-missing bootstrap = total obs. (_FREQ_) - nmiss ;
	do i=1 to dim(varlist);
		varlist(i) = _FREQ_ - (varlist(i));
	end;
	run;
	
	proc transpose data=_nboot out=_nboott(rename=(COL1=nboot)) prefix=COL ;
	by p ;
	var t: ;
	run;
	 
/* calculate SD of sampling distribution (bootstrap SE)  */
	proc means data=survivalplotbs std noprint;
	by p ;
	var t: ;
	output out=_std(drop=_TYPE_) std=;
	run;

/* transpose both estimate and boostrap SE */
	proc transpose data=_std out=_stdt(rename=(COL1=se)) prefix=COL ;
	by p ;
	var t: ;
	run;
	 
	proc transpose data=survivalplot0 out=_estimate (rename=(estimate1=estimate)) prefix=estimate;
	by p ;
	var t: ;
	run;
 
/* merge estimate and boostrap SE, calculate 95ci */
	data survivalplotf (rename=(_NAME_=name));
	merge _estimate _stdt _nboott;
	by p ;
	* calculate 95ci, pvalue;
	alpha = 1 - (95/100);
	one_minus_half_alpha = 1 - alpha/2;
	t_quant = quantile('t', one_minus_half_alpha, nboot-1);
	lcl95 = estimate - t_quant * se;
	ucl95 = estimate + t_quant * se;
	if (se ne 0) then do;
		tvalue =abs( ( estimate - 1 ) / se );
		format pvalue PVALUE6.4;
		pvalue = 2 * (1 - probt(tvalue, nboot-1 ) );
	end;
	label pvalue = "2-sided pvalue (null H value=1)"
		  _NAME_= " "
		  nboot = " ";
	run;

/* save data */
	data reslib.survivalplotf(rename=(hefi2019=HEFI2019_TOTAL_SCORE));
	retain p hefi2019 name time_to_primary ;
		merge survivalplotf reslib.Inexposure(keep=p hefi2019);
			by p;
	* add time to primary outcome as numerical value ;
	time_to_primary=input(compress(name,,'a'),10.);
	run;

 /*************************************************************************/
 /*                                                                       */
 /* SURVIVAL CONTRASTS: BOOTSTRAP VARIANCE ESTIMATION                     */
 /*                                                                       */
 /*************************************************************************/

/* separate original and bootstrap replicates */
	data survivalcontrast0 survivalcontrastbs;
		set reslib.survivalplot_w2;
	* scale Pr(outcome) and risk difference to percentage point;
	array values(*) p: rd: ;
	do i=1 to dim(values);
		values(i) = values(i)*100;
	end;
	if replicate=0 then output survivalcontrast0;
		else output survivalcontrastbs;
	run;

 /************************************************/
 /*       Pr(outcome) and risk difference        */
 /************************************************/

	proc sort data=survivalcontrastbs;
		by time_to_primary;
	run;
	 
	proc sort data=survivalcontrast0;
		by time_to_primary;
	run;

/* output number of successful (non-missing) bootstraps */
	proc means data=survivalcontrastbs nmiss noprint;
	by time_to_primary ;
	var p: rd: ;
	output out=_nboot(drop=_TYPE_) nmiss=;
	run;
	 
	data _nboot;
	set _nboot;
	array varlist (*) p: rd: ;
	* Number of non-missing bootstrap = total obs. (_FREQ_) - nmiss ;
	do i=1 to dim(varlist);
	varlist(i) = _FREQ_ - (varlist(i));
	end;
	run;

	proc transpose data=_nboot out=_nboott(rename=(COL1=nboot)) prefix=COL ;
	by time_to_primary ;
	var p: rd: ;
	run;

/* calculate SD of sampling distribution (bootstrap SE)  */
	proc means data=survivalcontrastbs std noprint;
	by time_to_primary ;
	var p: rd: ;
	output out=_std(drop=_TYPE_) std=;
	run;

/* transpose both estimate and boostrap SE */
	proc transpose data=_std out=_stdt(rename=(COL1=se)) prefix=COL ;
	by time_to_primary ;
	var p: rd: ;
	run;
 
	proc transpose data=survivalcontrast0 out=_estimate (rename=(estimate1=estimate)) prefix=estimate;
	by time_to_primary ;
	var p: rd: ;
	run;
	
/* merge estimate and boostrap SE, calculate 95ci */
	data survivalcontrastf1 (rename=(_NAME_=name));
	merge _estimate _stdt _nboott;
	by time_to_primary ;
	
	* calculate 95ci, pvalue;
	alpha = 1 - (95/100);
	one_minus_half_alpha = 1 - alpha/2;
	t_quant = quantile('t', one_minus_half_alpha, nboot-2);
	lcl95 = estimate - t_quant * se;
	ucl95 = estimate + t_quant * se;
	if (se ne 0) then do;
		tvalue =abs( ( estimate - 0 ) / se );
		format pvalue PVALUE6.4;
		pvalue = 2 * (1 - probt(tvalue, nboot-2 ) );
	end;
	
	label pvalue = "2-sided pvalue (null H value=0)" _NAME_= " " nboot = " ";
	run;

 /************************************************/
 /*         Risk ratio (relative risks)          */
 /************************************************/

/* output number of successful (non-missing) bootstraps */
	proc means data=survivalcontrastbs nmiss noprint;
	by time_to_primary ;
	var rr: ;
	output out=_nboot(drop=_TYPE_) nmiss=;
	run;
	 
	data _nboot;
	set _nboot;
	array varlist (*) rr: ;
	* Number of non-missing bootstrap = total obs. (_FREQ_) - nmiss ;
	do i=1 to dim(varlist);
	varlist(i) = _FREQ_ - (varlist(i));
	end;
	run;
 
	proc transpose data=_nboot out=_nboott(rename=(COL1=nboot)) prefix=COL ;
	by time_to_primary ;
	var rr: ;
	run;

/* calculate SD of sampling distribution (bootstrap SE)  */
	proc means data=survivalcontrastbs std noprint;
	by time_to_primary ;
	var rr: ;
	output out=_std(drop=_TYPE_) std=;
	run;

/* transpose both estimate and boostrap SE */
	proc transpose data=_std out=_stdt(rename=(COL1=se)) prefix=COL ;
	by time_to_primary ;
	var rr: ;
	run;
	 
	proc transpose data=survivalcontrast0 out=_estimate (rename=(estimate1=estimate)) prefix=estimate;
	by time_to_primary ;
	var rr: ;
	run;

/* merge estimate and boostrap SE, calculate 95ci */
	data survivalcontrastf2 (rename=(_NAME_=name));
	merge _estimate _stdt _nboott;
	by time_to_primary ;
	* calculate 95ci, pvalue;
	alpha = 1 - (95/100);
	one_minus_half_alpha = 1 - alpha/2;
	t_quant = quantile('t', one_minus_half_alpha, nboot-2);
	lcl95 = estimate - t_quant * se;
	ucl95 = estimate + t_quant * se;
	if (se ne 0) then do;
		tvalue =abs( ( estimate - 1 ) / se );
		format pvalue PVALUE6.4;
		pvalue = 2 * (1 - probt(tvalue, nboot-2 ) );
	end;
	label pvalue = "2-sided pvalue (null H value=1)" _NAME_= " " nboot = " ";
	run;
	 
 /*************************************************************************/
 /*                                                                       */
 /* SURVIVAL CONTRASTS: APPEND AND SAVE DATA                              */
 /*                                                                       */
 /*************************************************************************/

	data reslib.survivalcontrastf;
	retain type time_to_primary name p;
		set survivalcontrastf1 survivalcontrastf2;
	length type $ 18;
	* classify estimates;
		if index(name,'RR') > 0 then type="Risk ratio";
		else if index(name,'RD') > 0 then type="Risk difference";
		else type= "P(outcome)";
	* add percentile values for p(outcome);
	if type="P(outcome)" then do;
		p = input(compress(name,,'a'),10.);
	end;
	else p =.;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*                              END OF CODE                              */
 /*                                                                       */
 /*************************************************************************/