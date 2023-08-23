# new parameter names august 21, 2023

##############
# plot rhat diagnostic
pdf(file = paste(wd,'rhat_',current.chain,'.pdf',sep=''))
mcmc_rhat(rhats)
graphics.off()

# plot selectivities
pdf(file = paste(wd,'selectivies_',current.chain,'.pdf',sep=''))
par(cex=1.3)
plot(names(krill_mpd$Fsh_lf_Sel$Values),krill_mpd$Fsh_lf_Sel$Values, 
  main='Estimate Krill Selectivity',xlab = 'Length (mm)',ylab = 'Length Selectivity',type='l',lwd=3,col='red')
  lines(names(krill_mpd$AMLR_lf_Sel$Values),krill_mpd$AMLR_lf_Sel$Values,col='blue',lwd=3,lty=1)
  lines(names(krill_mpd$AMLR_ac_Sel$Values),krill_mpd$AMLR_ac_Sel$Values,col='blue',lwd=3,lty=2)
  lines(names(krill_mpd$Fsh_ac_Sel$Values),krill_mpd$Fsh_ac_Sel$Values,col='red',lwd=3,lty=2)
  legend(6,0.4,lty=c(1,1,2,2),col=c('red','blue','blue','red'),
    c('Fishery_lfs','AMLR_lfs','AMLRacoustics','FisheryAcoustics'),lwd=3,cex=0.8)
graphics.off()

##############
# plot selectivities
pdf(file = paste(wd,'selectivies_',current.chain,'.pdf',sep=''))
par(cex=1.3)
plot(names(krill_mpd$Fsh_lf_Sel$Values),krill_mpd$Fsh_lf_Sel$Values, 
  main='Estimate Krill Selectivity',xlab = 'Length (mm)',ylab = 'Length Selectivity',type='l',lwd=3,col='red')
  lines(names(krill_mpd$AMLR_lf_Sel$Values),krill_mpd$AMLR_lf_Sel$Values,col='blue',lwd=3,lty=1)
  lines(names(krill_mpd$AMLR_ac_Sel$Values),krill_mpd$AMLR_ac_Sel$Values,col='blue',lwd=3,lty=2)
  lines(names(krill_mpd$Fsh_ac_Sel$Values),krill_mpd$Fsh_ac_Sel$Values,col='red',lwd=3,lty=2)
  legend(6,0.4,lty=c(1,1,2,2),col=c('red','blue','blue','red'),
    c('Fishery_lfs','AMLR_lfs','AMLRacoustics','FisheryAcoustics'),lwd=3,cex=0.8)
graphics.off()

##############
# plot median and 95% CIs for mcmc estimates of ssbs, overlay mpd estimates of medians
pdf(file=paste(wd,'SSB_',current.chain,'.pdf',sep=''))
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

##############
# plot median and 95% CIs for mcmc estimates of recruitment
pdf(file=paste(wd,'Recruits_',current.chain,'.pdf',sep=''))
par(cex=1.2)
y.lim <-c(0,max(1,quantile_recruits_df[,3]))
yrs <- as.numeric(c(substr(rownames(quantile_recruits_df),10,13)))
plot(yrs,quantile_recruits_df[,2],type='l',lwd=3,ylim=y.lim,
  xlab= 'Year',ylab='Recruits')
lines(yrs,quantile_recruits_df[,1],lwd=3)
lines(yrs,quantile_recruits_df[,3],lwd=3)
polygon(c(yrs,rev(yrs)),
  c(quantile_recruits_df[,2],rev(quantile_recruits_df[,3])),col='light grey')
polygon(c(yrs,rev(yrs)),
  c(quantile_recruits_df[,2],rev(quantile_recruits_df[,1])),col='light grey')
lines(yrs,quantile_recruits_df[,2],lwd=3)
#abline(h=median(cas2_tab$summary$values[[1]]),lwd=2,lty=2)
points(yrs,krill_mpd$Recruitment$recruits,col='red',pch=19)
#segments(
#  yrs,quantile_recruit_mpd_df[,1], 
#  yrs,,quantile_recruit_mpd_df[,3], 
#  lwd=2,col='red')
graphics.off()
