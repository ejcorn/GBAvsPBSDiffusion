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

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months

Grp.mean <- colMeans(Mice,na.rm = T)

# Fit time scaling parameter on average of all mice
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- log(Grp.mean,base=10)
seed.fits <- matrix(NA,n.regions) # preallocate vector to store fits from each seed region

for(S in 1:n.regions){
  seed.name <- region.names[S]
  print(paste('Seed:',seed.name))
  Xo <- get.Xo(region.names,seed.name) # seed pathology in each region one at a time
  dat <- c.fit.r(log.path,L.out,Xo,c.rng) # fit time constant
  seed.fits[S] <- dat$r[1] # store correlations from month 1
}

rownames(seed.fits) <- region.names
save(seed.fits, file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))