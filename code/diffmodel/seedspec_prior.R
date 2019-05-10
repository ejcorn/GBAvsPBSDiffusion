# This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel_priorfit/seedspec/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')

#################
### Load data ###
#################

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression
Synuclein <- norm.synuclein(Synuclein,1)

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months
tp <- c(1,3,6)
Mice.prior <- lapply(tp[2:3], function(M) prior.data[prior.data$MPI == M,-1])

# join pathology measurements
Mice <- c(list(Mice),Mice.prior)
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

# Fit time scaling parameter on average of all mice
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)
L.out <- get.Lout(W,as.numeric(Synuclein)) # compute out-degreee Laplacian weighted by Snca

c.rng <- seq(params$c.min,2*params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
seed.fits <- matrix(NA,n.regions) # preallocate vector to store fits from each seed region

for(S in 1:n.regions){
  seed.name <- region.names[S]
  print(paste('Seed:',seed.name))
  Xo <- get.Xo(region.names,seed.name) # seed pathology in each region one at a time
  dat <- c.fit.prior.r(log.path,L.out,Xo,tp,c.rng) # fit time constant
  seed.fits[S] <- dat$r[1] # store correlations from month 1
}

rownames(seed.fits) <- region.names
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

ggsave(p.null.seeds,filename=paste(savedir,grp,'iCPSeedSpecificityPriorFit.pdf',sep=''),width=2,height=3,units='in')