# Casal2-krill-model-mcmc
MCMC sampling of recruitment, spawning stock biomass and the application of CCAMLR decision rules are added to the Casal2 model for Antarctic krill in Subarea 48.1
 
This repository updates the pilot Casal2 krill model described in 'https://github.com/us-amlr/Casal2-krill-model'. The update uses 'r4Casal2_krill.r' to evaluate derived parameters for spawning stock biomass and calculates the CCAMLR 'depletion' and 'escapement' decision rules in order to define precautionary yield for the Antarctic krill fishery in CCAMLR Subarea 48.1. After 'r4Casal2_krill.r' is applied to the Casal2 output files the results can be plotted using 'r4Casal2_krilPLOTS.r'.

The repository contains four alternative configurations of the krill model, 'combined non-simplex', 'combined simplex', 'separate non-simplex' and 'separate simplex'. These employ two different methods of estimating annual recruitment multipliers - simplex transformation and non-simplex - and two alternative approaches to modeling future projected catches, either with all model years combined into a single ‘@process Instantaneous_Mortality’ block, or with future catch years separated from past catch years into a second ‘@project future_catch’ block.

To run the configurationss I place the 'config' directory in one of the four github configuration directories into a local 'biom' directory. The 'biom' directory also gets the unzipped 'Casal.exe.zip' and 'dll.zip' files. The 'biom' directory will also contain the output files after a model run. To run multiple MCMC chains (I run three per configuration) put the same 'config', 'Casal2.exe' and '.dll' files in differently labeled directories such as 'a', 'b' and 'c'. I label the directories with the projected future catches as well, i.e. '0K.a' for the first directory with 0 tonnes of future catches. Each directory with identical input files will start from a different random seed, representing a separate mcmc chain. The default settings will take about 25 minutes to run for each chain. Model diagnostics such as the 'rhat'  included with 'r4Casal2_krill.r' and those in the R package 'coda' will improve with longer iterations of MCMC samples. 

Run the 3 chains with all recruitment_multipliers in population.csl2 equal to 1.0 (each model will start at a different random seed) from the Windows console using the following commands:

cd 'localpath'\biom

casal2 -r >> krill.r.txt

casal2 -e >> krill.e.txt

casal2 -m >> krill_mcmc.m.txt

casal2 -r -i samples.1 -t >> krill_samples.txt

where 'localpath' is the local path to each of the three directories containing 'biom' subdirectories with the Casal2 files.


The runs with all recruitment_multipliers equal to 1 will produce a constant number of recruits and spawning stock biomass in the future projections. This will not provide a test of the CCAMLR depletion rule because there will be no variability in the future projections. In order to assess the depletion rule, recruitment fluctuations estimated from past data can be applied to the twenty year projection in a new model by starting the 20-year future projection with the recruitment variability estimates from the AMLR summer survey data from 1992-2011. Fishery catches during this period were generally small, under 150,000 tonnes/year.

In a second set of runs representing the population in future years with catches, enter the recruitment variability from the first 0 catch year for 1992-2011 ('recruitment_multipliers[17:36]') as the 'recruitment_multipliers' for the 20 projection years in 'population.csl2'. 

Once the runs are finished the mcmc samples of the spawning biomass posterior distribution can be sampled and plotted, and the CCAMLR Decision rules applied, using 'r4Casal2_krill.r'. Paths in 'r4Casal2_krill.r' will need to be set to the local file structure. MPD values may be plotted using 'https://github.com/us-amlr/Casal2-krill-model/blob/main/R_krill_figs.rmd'.

The rhats directory illustrates the effect of increasing MCMC sample size from 110,000 (25 minutes run time) to 3,300,000 (17 hours run time) on the Rhats diagnostic plot from the rstan library.
