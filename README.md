# Casal2-krill-model-mcmc
MCMC analysis of spawning stock biomass uncertainty and the application of CCAMLR decision rules are added to the Casal2 model for Antarctic krill in Subarea 48.1
 
This repository updates the pilot Casal2 krill model described in 'https://github.com/us-amlr/Casal2-krill-model'. The update uses 'r4Casal2_krill.r' to evaluate derived parameters for spawning stock biomass and calculates the CCAMLR 'depletion' and 'escapement' decision rules in order to define precautionary yield for the Antarctic krill fishery in CCAMLR Subarea 48.1.

'SSB_0K.a.pdf' is a plot of the spawning stock biomass estimates derived from MCMC sampling during the period from 1976-2021 with reported fishery catches in the model. All recruitment multipiers were initialized to 1.0 and zero future catches were projected forward from 2022 to 2041. 'SSB_600K.a.pdf' shows the model spawning stock biomass with 600 tonnes annual yields projected forward from 2022 to 2041. Recruitment multipliers for the projection years were taken from the estimates for the observed years 1992 to 2011 (SSB_600K.a.pdf). These model outputs can be recreated using the files in this repository.

The repository contains two alternative versions of the krill model. The output files for the first model are in the 'biom' directory and the input files are in 'biom/config'. The second model alternative uses the simplex transformation to estimate recruitment multipliers. Inputs and outputs for this configuration are found in the 'simplex_rec_multipliers' directory.

The two configurations produce different results. Because the data and data weightings are the same for the two configurations, they can be compared using AIC. AIC = 2K -2ln(L) where K is the number of estimated parammeters and ln(L) is the log-likelihood. 

In each of the two configurations the negative log-likelihood is reported in the 'krill.e.txt' file labeled as 'total_score'. This score is -181.115 for the non-simplex model and -87.5165 for the simplex model. The non-simplex model estimates 55 parameters and the non-simplex model estimates 53. Hence the non-simplex AIC score is:

(2*55) -(2*181.115) = -252.23 

and the simplex score is:

(2*53) -(2*87.5165) = -69.033.

This difference of over 180 AIC between the configurations indicates the non-simplex model is a much better representation of these data.

To run the models place three sets of the files in the 'biom' directory in differently labeled directories such as 'a', 'b' and 'c'. I label the directories with the projected future catches as well, i.e. '0K.a' for the first directory with 0 tonnes of future catches.

Put 'casal2.exe', 'config.csl2' and the four compressed files in the 'dlls' directory in the three directories containing the 'biom' subdirectories.
Run the 3 chains with all recruitment_multipliers in population.csl2 equal to 1.0 (each model will start at a different random seed) from the console using the following commands:

cd 'localpath'\biom

casal2 -r >> krill.r.txt

casal2 -e >> krill.e.txt

casal2 -m >> krill_mcmc.m.txt

casal2 -r -i samples.1 -t >> krill_samples.txt

where 'localpath' is the local path to each of the directories containing three 'biom' subdirectories and other files.

Each directory with identical input files will start from a different random seed, representing a separate mcmc chain. The default settings will take about 25 minutes to run for each chain. Model diagnostics such as the 'rhat'  included with 'r4Casal2_krill.r' and those in the R package 'coda' will improve with longer iterations. 

The runs with all recruitment_multipliers equal to 1 will produce a constant number of recruits and spawning stock biomass in the future projections. This will not provide a test of the CCAMLR depletion rule because there will be no variability in the future projections. In order to assess the depletion rule, recruitment fluctuations estimated from past data can be applied to the twenty year projection in a new model by starting the 20-year future projection with the recruitment variability estimates from the AMLR summer survey data from 1992-2011. Fishery catches during this period were generally small, under 150,000 tonnes/year.

In a second set of runs, enter the recruitment variability from the first 0 catch year for 1992-2011 ('recruitment_multipliers[17:36]') as the 'recruitment_multipliers' for the last 20 projection years in 'population.csl2'. 

Once the runs are finished the mcmc samples of the spawning biomass posterior distribution can be sampled and plotted, and the CCAMLR Decision rules applied, using 'r4Casal2_krill.r'. Paths in 'r4Casal2_krill.r' will need to be set to the local file structure. MPD values may be plotted using 'https://github.com/us-amlr/Casal2-krill-model/blob/main/R_krill_figs.rmd'.

The rhats directory illustrates the effect of increasing MCMC sample size from 110,000 (25 minutes run time) to 11,000,000 (47 hours run time) on the Rhats diagnostic plot from the rstan library and on the CCAMLR decision rule for escapement.
