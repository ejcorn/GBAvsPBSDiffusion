#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'PBSvsCBE/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

c.rng <- seq(0.1,60,length.out = 2000) # scaling parameter range
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu

predicted.path.c <- lapply(c.rng, function(c) predict.Lout(L.out,Xo,c,t=1)) # predictions for a bunch of c-values
save(c.rng,predicted.path.c,file=paste(savedir,'XtForCRng.RData',sep=''))