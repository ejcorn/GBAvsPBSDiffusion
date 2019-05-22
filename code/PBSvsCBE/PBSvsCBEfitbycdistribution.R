# this script will compute model fits for each mouse at multiple time constant values
# then plot the variance of the fits at each value of c

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'PBSvsCBE/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

# get path data for group
Mice <- path.data[path.data$Condition == grp,-1]

# fits from iCPU
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
n.mice <- nrow(Mice)
tf <- 0.5 # training fraction
c.train.Grp <- r.Grp <- rep(NA,n.mice)

load(file=paste(savedir,'XtForCRng.RData',sep=''))

fit.list <- list()
for(MOUSE in 1:n.mice){
	print(paste('Mouse',MOUSE))
	# get path data for each mouse
	path.data.train.Grp <- colMeans(Mice[MOUSE,])
	log.path.train <- log(path.data.train.Grp,base=10)
	mask <- path.data.train.Grp != 0 # identify regions with no path --> log = -Inf
	# compute fit between each mouse's data and predicted path for range of c-values
	fit.by.c <- sapply(predicted.path.c, function(X) cor(log(X,base=10)[mask],log.path.train[mask]))
	fit.list[[MOUSE]] <- fit.by.c
	# store c-value with best fit
	c.train.Grp[MOUSE] <- c.rng[which.max(fit.by.c)]
	# store fit at that c-value
	r.Grp[MOUSE] <- max(fit.by.c)
}

X <- do.call('rbind',fit.list) # compile fits into Mouse-by-c matrix
X.sd <- colSD(X)
X.mean <- colMeans(X)
df <- data.frame(m=X.mean,s=X.sd,x=c.rng)
df$grp <- rep(grp,nrow(df))

save(df,file=paste(savedir,grp,'FitsByCByMouse.RData',sep=''))