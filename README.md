# Casal2-krill-model-update
 Adds analysis of spawning stock biomass and CCAMLR decision rules to 'us-amlr/Casal2-krill-model'
 
This repository contains updates to the Casal2 krill model described in 'https://github.com/us-amlr/Casal2-krill-model'. These updates use 'r4Casal2_krill.r' to evaluate derived parameters for spawning stock biomass and calculate the CCAMLR 'depletion' and 'escapement' decision rules in order to define precautionary yield.

To run the model place three sets of the files in the 'biom' directory in differently labeled directories such as 'a', 'b' and 'c'. I label the directories with the catch as well, i.e. '0K.a' for the first directory with 0 tonnes catch.
Put 'casal2.exe', 'config.csl2' and the four .dll files in the three directories containing the 'biom' subdirectory.
Run the 3 models with all recruitment_multipliers in population.csl2 equal to 1.0 (each model will start at a different random seed) from the console using the following commands:

cd 'localpath'\biom
casal2 -r >> krill.r.txt
casal2 -e >> krill.e.txt
casal2 -m >> krill_mcmc.m.txt 
casal2 -r -i samples.1 -t >> krill_samples.txt

where 'localpath' is the local path to each of the directories containing three 'biom' subdirectories and other files.
The default mcmc settings will take about 25 minutes to run. Model diagnostics such as the 'rhat'  included with 'r4Casal2_krill.r' and those in the R package 'coda' will improve with longer iterations. 

The runs with all recruitment_multipliers equal to 1 will produce a constant number of recruits and spawning stock biomass in the future projections. This will not provide a test of the CCAMLR depletion rule because there will be no variability in the future projections. In order to assess the depletion rule, recruitment fluctuations estimated from past data can be applied to the twenty year projection in a new model by starting the 20-year future projection with the recruitment variability estimates from the AMLR summer survey data from 1992-2011. Fishery catches during this period were generally small, under 150,000 tonnes/year.

In a second set of runs, enter the recruitment variability from the first 0 catch year for 1992-2011 ('recruitment_multipliers[17:36]') as the 'recruitment_multipliers' for the last 20 projection years in 'population.csl2'.

Once the runs are finished the mcmc samples of the spawning biomass posterior distribution can be sampled and plotted, and the CCAMLR Decision rules applied, using 'r4Casal2_krill.r'. MPD values may be plotted using 'https://github.com/us-amlr/Casal2-krill-modelR_krill_figs.rmd'.
