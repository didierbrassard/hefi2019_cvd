Adherence to Canada’s Food Guide and cardiovascular disease
================

This repository presents the analysis code for the study *Greater
adherence to the 2019 Canada’s Food Guide recommendations on healthy
food choices reduces the risk of cardiovascular disease in adults: a
prospective analysis of UK Biobank data* published in the American
Journal of Clinical Nutrition (Brassard et al., 2022a).

The general objective was to estimate the effect of adherence to
Canada’s Food Guide recommendations on healthy food choices, measured
using the HEFI-2019 (Brassard et al, 2022b, 2022c), on cardiovascular
disease (CVD) risk in the UK Biobank prospective study.

# 1. Requirements and file structure

## 1.1 Data and software

The baseline covariate data (November 2020), original 24-h recall
dietary intake data (July 2021), and outcome data (April 2021) of the UK
Biobank were used, under application \#25205.

The main analyses were executed in SAS (v9.4; maintenance release
9.04.01M7P08052020) in a Windows 10 64-bits environment. The manuscript
file was successfully executed using RMarkdown and R (version 4.1.3) on
macOS Big Sur 11.6.8 (64-bit).

## 1.2 SAS macros

SAS macros made by the National Cancer Institute (NCI) were used to
perform measurement error correction (Zhang et al., 2011). The macros
are [available on the NCI
website](https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0),
but are included in the present repository as well in the `/Macros/`
folder.

The HEFI-2019 scoring algorithm SAS macro was used and it is also
[available online](https://github.com/didierbrassard/hefi2019). Of note,
the version in the present repository includes a very minor
modification. Since beverages were reported in volume (ml) in the UK
Biobank, the amount of unsweetened protein-rich beverages for 1
reference amount was set at 250 ml in the present study instead of 258g
in the original macro (`%let probev_gram_per_RA = 250;` L143 of the
modified macro).

## 1.3 Structure

Key steps of the main analysis are presented in separate `.sas` files.
Each file would need to be executed in sequential order for successful
execution. Raw UK Biobank file are not available publicly. However, the
UK Biobank is an open access resource. Data described in the manuscript
and code book are available upon request by registering and applying
online (<http://www.ukbiobank.ac.uk/register-apply>).

# 2. Brief description

Complete details about analysis are provided in the article as well as
in the online supplementary material.

**[Step 1. Distribution of HEFI-2019
score](01_HEFI19CVD_Usual_intake_distribution.sas)**: estimation of
population-level usual dietary intakes with the NCI multivariate methods
to obtain distribution of HEFI-2019 score.

**[Step 2. Participant-level
simulation](02_HEFI19CVD_Usual_intake_simulation.sas)**: simulation of
usual dietary intakes at the participant level with the NCI multivariate
method to obtain measurement-error-corrected regression coefficients in
outcome models.

**[Step 3. Inverse probability
weighting](03_HEFI19CVD_Inverse_probability_weighting.sas)**: estimation
of inverse probability of “treatment” weights, inverse probability of
censoring weights, and combined weights to adjust for confounding and
censoring.

**[Step 4. Outcome model](04_HEFI19CVD_Outcome_model.sas)**: estimation
of survival curves using energy-adjusted relationship between the total
HEFI-2019 score based on usual intakes and incident CVD with Cox
proportional hazards regression models.

**[Step 5. Bootstrap variance
estimation](05_HEFI19CVD_Bootstrap_variance.sas)**: parametric bootstrap
variance estimation using the 250 resampling to obtain confidence
intervals.

# 3. References

Brassard D, Manikpurage Hasanga D, Thériault S, et al. Greater adherence
to the 2019 Canada’s Food Guide recommendations on healthy food choices
reduces the risk of cardiovascular disease in adults: a prospective
analysis of UK Biobank data. Am J Clin Nutr 2022a. doi: nqac256

Brassard D, Elvidge Munene LA, St-Pierre S, et al. Development of the
Healthy Eating Food Index (HEFI)-2019 measuring adherence to Canada’s
Food Guide 2019 recommendations on healthy food choices. Appl Physiol
Nutr Metab 2022b;47:595-610. doi:
[10.1139/apnm-2021-0415](https://doi.org/10.1139/apnm-2021-0415).

Brassard D, Elvidge Munene LA, St-Pierre S, et al. Evaluation of the
Healthy Eating Food Index (HEFI)-2019 measuring adherence to Canada’s
Food Guide 2019 recommendations on healthy food choices. Appl Physiol
Nutr Metab 2022c;47(5):582-94. doi:
[10.1139/apnm-2021-0416](https://doi.org/10.1139/apnm-2021-0416).

Zhang S, Midthune D, Guenther PM, et al. A New Multivariate Measurement
Error Model with Zero-Inflated Dietary Data, and Its Application to
Dietary Assessment. Ann Appl Stat 2011;5(2B):1456-87. doi:
[10.1214/10-AOAS446](https://doi.org/10.1214/10-AOAS446).
