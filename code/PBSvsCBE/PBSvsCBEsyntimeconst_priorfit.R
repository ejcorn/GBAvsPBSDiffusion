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

# get data from all mice for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months
tp <- c(1,3,6)
Mice.prior <- lapply(tp[2:3], function(M) prior.data[prior.data$MPI == M,-1])
Grp.mean.prior <- lapply(Mice.prior, function(x) matrix(colMeans(x,na.rm = T),nrow=1))

# join pathology measurements
Mice <- c(list(Mice),Grp.mean.prior)

# Fit time scaling parameter on average of all mice
c.rng <- seq(0.5,1.5,length.out = params$c.n*0.75) # scaling parameter
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu

n.reps <- 100
tf <- 'SRS' # simple random sample with replacement
c.train.Grp <- rep(NA,n.reps)
r.Grp <- list()

for(REP in 1:n.reps){
	print(paste('Rep',REP))
	# generate indices for bootstrapping
	train.idx.Grp <- list(sample(1:nrow(Mice[[1]]),replace=T),1,1)
	# split into non-overlapping training and testing sets
	path.data.train.Grp <- mapply(function(X,t.idx) {list(as.matrix(X[t.idx,]))}, X=Mice,t.idx=train.idx.Grp)
	# compute train and test means
	path.data.train.Grp[[1]] <- colMeans(path.data.train.Grp[[1]],na.rm=T)
	# compute log
	log.path.train <- lapply(path.data.train.Grp, function(x) log(x,base=10))
	# fit time scale parameter on training data
	list[c.train.Grp[REP],r.Grp[[REP]]] <- c.fit.prior.r(log.path.train,L.out,Xo,tp,c.rng)
}

save(c.train.Grp,r.Grp,file = paste(savedir,grp,'Syntimeconstants_priorfitTF',tf,'.RData',sep=''))