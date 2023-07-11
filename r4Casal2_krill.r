# modified from C. Marsh code at https://niwafisheriesmodelling.github.io/r4Casal2/mcmc.html)
# cloned and modified r4Casal2 from 'https://github.com/NIWAFisheriesModelling/r4Casal2'

 rm(list=ls())
 lib.path <- 'c:/zot/Casal2/2023/'
 path <- paste(lib.path,'R-Libraries/R/',sep='')
 r.scrpts <- dir(path)
 for(i in 1:length(r.scrpts))
   source(paste(path,r.scrpts[i],sep=''))
library(knitr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(r4Casal2)

ch.names <- c('500K_long.a','500K_long.b','500K_long.c')           # identically configured model names for mcmc chains
wd <- 'C:/zot/Casal2/2023/7jul/jul10/'              # working directory for mcmc chains
current.chain <- '500K_long.c'                           # chain to be plotted
wd2 <- paste(wd,current.chain,'/biom/',sep='')      # working directory for Casal2 input files

setwd(wd)

mcmc.nm <- mcmc_post <- list()
for(i.ch in 1:(length(ch.names))){
  mcmc.nm[[i.ch]] = extract.mcmc(path = paste(wd,ch.names[i.ch],'/biom',sep=''), 
    samples.file = "samples.1", objectives.file = "objectives.1")
  mcmc.nm[[i.ch]]$chain = as.character(i.ch)
  mcmc_post[[i.ch]] = mcmc.nm[[i.ch]] %>% filter(state == "mcmc")
  }

mcmc_all <- rbind(mcmc_post[[1]],mcmc_post[[2]],mcmc_post[[3]])
mcmc_non_burn_in = mcmc_all %>% filter(state == "mcmc")
n_posterior_samples =  nrow(mcmc_post[[1]])+nrow(mcmc_post[[2]])+nrow(mcmc_post[[3]])
pars = colnames(mcmc_post[[1]][,12:(ncol(mcmc_non_burn_in) - 1)])
iters = max(nrow(mcmc_post[[1]]),nrow(mcmc_post[[2]]),nrow(mcmc_post[[3]]))
bayes_array = array(dim = c(iters, 3, length(pars)), dimnames = list(1:iters, 1:3, pars))
bayes_array[1:nrow(mcmc_post[[1]]),1,] = as.matrix(mcmc_post[[1]][,12:(ncol(mcmc_non_burn_in) - 1)])
bayes_array[1:nrow(mcmc_post[[2]]),2,] = as.matrix(mcmc_post[[2]][,12:(ncol(mcmc_non_burn_in) - 1)])
bayes_array[1:nrow(mcmc_post[[3]]),3,] = as.matrix(mcmc_post[[3]][,12:(ncol(mcmc_non_burn_in) - 1)])
min_cutoff = min(nrow(mcmc_post[[1]]),nrow(mcmc_post[[2]]),nrow(mcmc_post[[3]]))
bayes_array = bayes_array[1:min_cutoff, ,]

library(purrr)
library(bayesplot)
# library(posterior)
library(rstan) # required for function Rhat

rhats = apply(bayes_array, MARGIN = 3, Rhat)
n_eff_bulk = apply(bayes_array, MARGIN = 3, ess_bulk)
n_eff_tail = apply(bayes_array, MARGIN = 3, ess_tail)

pdf(file = paste(wd,'mcmc_rhat_',current.chain,'.pdf',sep=''))
mcmc_rhat(rhats)
graphics.off()

# set probabilities to calculate from mcmc samples
p <- c(0.025, 0.5, 0.975) ## confidence intervals
p_names <- map_chr(p, ~paste0(.x*100, "%"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  rlang::set_names(nm = c("low", "mid", "upp"))

# extract derived SSB values from tabular mcmc estimates
#cas2_file_name = system.file("extdata", "tabular.log", package = "r4Casal2", mustWork = TRUE) # original example
cas2_file_name = paste(wd2,'krill_samples.txt',sep='') # this is samples.1 produced from casal2 -m,
                                                       # casal2 -r -i samples.1 -t >> krill_samples.txt
cas2_tab = extract.tabular(file = cas2_file_name, quiet = T)
## cut off burn-in the first 50 samples
#cas2_tab = burn.in.tabular(cas2_tab, Row = 50)

# get selectivities
selectivity_df = get_selectivities(cas2_tab)
selectivity_df$selectivity_label <- selectivity_df$report_label # dhk hack to get dplyr plotting code to work
quantile_selectivity_df = selectivity_df %>% 
  group_by(bin, selectivity_label) %>% 
  summarize_at(vars(selectivity), p_funs)


# get annual spawning stock biomasses
ssbs = get_derived_quanitites(model = cas2_tab)
ssbs$years[ssbs$years == "initialisation_phase_1"] = 1976
#ssbs$years <- as.numeric(ssbs$years)
ssbs <- cbind(ssbs,B0=rep(cas2_tab$summary$values[[1]],67)) # add B0 column to ssbs

# calculate mcmc confidence intervals for spawning biomasses
quantile_ssb_df = ssbs %>% 
  group_by(years, dq_label) %>% 
  summarize_at(vars(values), p_funs)
quantile_ssb_df$years = as.numeric(quantile_ssb_df$years)

# calculate mpd estimates of median spawning biomasses
# (note that get_derived_quanitites does not produce uncertainty estimates for mpds)
krill_mpd = extract.mpd(paste(wd2,'krill.e.txt',sep='')) # extract mpd values at this location
ssbs_mpd = get_derived_quanitites(model = krill_mpd)
quantile_ssb_mpd_df = ssbs_mpd %>% 
  group_by(years, dq_label) %>% 
  summarize_at(vars(values), p_funs)
quantile_ssb_mpd_df$years = as.numeric(quantile_ssb_mpd_df$years)

# plot selectivities
pdf(file = paste(wd,'selectivies_',current.chain,'.pdf',sep=''))
par(cex=1.3)
plot(names(krill_mpd$krillFSel$Values),krill_mpd$krillFSel$Values, 
  main='Estimate Krill Selectivity',xlab = 'Length (mm)',ylab = 'Length Selectivity',type='l',lwd=3,col='red')
  lines(names(krill_mpd$trawlSel$Values),krill_mpd$trawlSel$Values,col='blue',lwd=3,lty=1)
  lines(names(krill_mpd$AMLR_Sel$Values),krill_mpd$AMLR_Sel$Values,col='blue',lwd=3,lty=2)
  lines(names(krill_mpd$FbiomSel$Values),krill_mpd$FbiomSel$Values,col='red',lwd=3,lty=2)
  legend(6,0.4,lty=c(1,1,2,2),col=c('red','blue','blue','red'),
    c('Fishery_lfs','AMLR_lfs','AMLRacoustics','FisheryAcoustics'),lwd=3,cex=0.8)
graphics.off()

# plot median and 95% CIs for mcmc estimates ssbs, overlay mpd estimates of medians
pdf(file=paste('SSB_',current.chain,'.pdf',sep=''))
par(cex=1.2)
y.lim <-c(0,max(1,quantile_ssb_df$upp))
plot(quantile_ssb_df$years,quantile_ssb_df$mid,type='l',lwd=3,ylim=y.lim,
  xlab= 'Year',ylab='SSB')
lines(quantile_ssb_df$years,quantile_ssb_df$low,lwd=3)
lines(quantile_ssb_df$years,quantile_ssb_df$upp,lwd=3)
polygon(c(quantile_ssb_df$years,rev(quantile_ssb_df$years)),
  c(quantile_ssb_df$mid,rev(quantile_ssb_df$upp)),col='light grey')
polygon(c(quantile_ssb_df$years,rev(quantile_ssb_df$years)),
  c(quantile_ssb_df$mid,rev(quantile_ssb_df$low)),col='light grey')
lines(quantile_ssb_df$years,quantile_ssb_df$mid,lwd=3)
abline(h=median(cas2_tab$summary$values[[1]]),lwd=2,lty=2)
points(quantile_ssb_mpd_df$years,quantile_ssb_mpd_df$mid,col='red',pch=19)
segments(
  quantile_ssb_mpd_df$years,quantile_ssb_mpd_df$low, 
  quantile_ssb_mpd_df$years,quantile_ssb_mpd_df$upp, 
  lwd=2,col='red')
graphics.off()

# Decision Rules
Gamma1 <- ssbs %>% group_by(years) %>% 
	  summarize(Dep=min(ssbs$values/ssbs$B0)) %>% 
	  summarize(Pr=mean(Dep < 0.2))
Gamma1

Gamma2 <- ssbs %>% 
	filter(years %in% max(years)) %>% 
	summarise(ssb=median(values),ssb0=median(cas2_tab$summary$values[[1]]))
Gamma2$Escapement<-Gamma2$ssb/Gamma2$ssb0
Gamma2

#The actual Gamma is the smallest of the two gammas:
GammaToUse<-which(c(Gamma1,Gamma2)==min(Gamma1,Gamma2)) #Which gamma is min?
if(length(GammaToUse)==2){GammaToUse=3} #when gamma1 and gamma2 are equal
OUT<-cbind(Gamma1,Gamma2,GammaToUse)
OUT

# Recruitment estimates
krill_mpd$Recruitment$standardised_recruitment_multipliers
krill_mpd$Recruitment$recruitment_multipliers
krill_mpd$Recruitment$true_ycs
krill_mpd$Recruitment$recruitment_multipliers[17:36]
