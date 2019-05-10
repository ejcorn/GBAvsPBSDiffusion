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
load(paste(params$opdir,'processed/Lout.RData',sep=''))

# get path data for group
Mice <- path.data[path.data$Condition == grp,-1]

# fits from iCPU
c.rng <- seq(params$c.min,params$c.max,length.out = 50) # scaling parameter range
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
n.reps <- 100
tf <- 0.5 # training fraction
c.train.Grp <- r.Grp <- rep(NA,n.reps)

for(REP in 1:n.reps){
	print(paste('Rep',REP))
	# mark training samples for each time point
	train.idx.Grp <- sample(1:nrow(Mice),size = tf*nrow(Mice),replace=FALSE)
	# compute mean for training sample
	path.data.train.Grp <- colMeans(Mice[train.idx.Grp,])
	#path.data.train.Grp <- colMeans(Mice[REP,])
	log.path.train <- log(path.data.train.Grp,base=10)	
	list[c.train.Grp[REP],r.Grp[REP],X] <- c.fit.r(log.path.train,L.out,Xo,c.rng)
}

save(c.train.Grp,r.Grp,file = paste(savedir,grp,'SyntimeconstantsTF',tf,'.RData',sep=''))