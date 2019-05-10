#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'PBSvsCBE/',sep='')
dir.create(savedir,recursive=T)
#hdir.create(paste(savedir,'roilevel',sep=''),recursive = T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

# get path data for group
Mice <- path.data[path.data$Condition == grp,-1]

# fits from iCPU
c.rng <- seq(params$c.min,params$c.max,length.out = 100) # scaling parameter range
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
n.mice <- nrow(Mice)
tf <- 0.5 # training fraction
c.train.Grp <- r.Grp <- rep(NA,n.mice)

for(MOUSE in 1:n.mice){
	print(paste('Mouse',MOUSE))
	# get path data for each mouse
	path.data.train.Grp <- colMeans(Mice[MOUSE,])
	log.path.train <- log(path.data.train.Grp,base=10)	
	list[c.train.Grp[MOUSE],r.Grp[MOUSE],X] <- c.fit.r(log.path.train,L.out,Xo,c.rng)
}

save(c.train.Grp,r.Grp,file = paste(savedir,grp,'TimeconstantsTF',tf,'.RData',sep=''))