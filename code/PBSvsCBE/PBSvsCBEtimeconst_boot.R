#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','predicted.path.c')))
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
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
n.reps <- 1000
tf <- 0.8 # training fraction
c.train.Grp <- r.Grp <- rep(NA,n.reps)

load(file=paste(savedir,'XtForCRng.RData',sep=''))

for(REP in 1:n.reps){
	print(paste('Rep',REP))
	# mark training samples for each time point
	train.idx.Grp <- sample(1:nrow(Mice),size = tf*nrow(Mice),replace=TRUE)
	# compute mean for training sample
	path.data.train.Grp <- colMeans(Mice[train.idx.Grp,])
	#path.data.train.Grp <- colMeans(Mice[REP,])
	log.path.train <- log(path.data.train.Grp,base=10)	
	mask <- path.data.train.Grp != 0 # identify regions with no path --> log = -Inf
	# compute fit between each sample's data and predicted path for range of c-values
	fit.by.c <- sapply(predicted.path.c, function(X) cor(log(X,base=10)[mask],log.path.train[mask]))
	# store c-value with best fit
	c.train.Grp[REP] <- c.rng[which.max(fit.by.c)]
	# store fit at that c-value
	r.Grp[REP] <- max(fit.by.c)
}

save(c.train.Grp,r.Grp,file = paste(savedir,grp,'SyntimeconstantsTF',tf,'.RData',sep=''))