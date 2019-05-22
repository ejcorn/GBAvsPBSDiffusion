# This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/seedspec/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')

#################
### Load data ###
#################

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))
# compute p-value as probability that other seeds predict observed path better than iCPu seed
null.cors <- seed.fits[region.names != 'iCP']
iCP.cor <- seed.fits[region.names == 'iCP']
seed.region.pval <- mean(iCP.cor < null.cors)
print(paste('Better fit than iCP:',region.names[which(seed.fits[region.names == 'iCP'] < seed.fits[region.names != 'iCP'])]))
month <- 'Month 1'
months <- rep(month,length(null.cors))
p.null.seeds <- ggplot() + 
  geom_jitter(aes(x=months,y = null.cors),color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1, position = position_dodge(width=1)) +
  geom_point(aes(x=month,y=as.numeric(iCP.cor)),shape = 18,color = 'black',size=2) + 
  geom_text(aes(x=month,y=0.8,label = paste('p =',signif(seed.region.pval,2))),size=2.5) +
  xlab('') + ylab('Fit') + ggtitle('Actual vs. Random Seed') +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8)) +
  theme(axis.text.x = element_text(size=8)) + theme(axis.text.y = element_text(size=8))
p.null.seeds

ggsave(p.null.seeds,filename=paste(savedir,grp,'iCPSeedSpecificity.pdf',sep=''),width=1.5,height=1.5,units='in')
