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

# load predicted values
X <- readMat(paste(savedir,'RegionalTrajectoriesCRange.mat',sep=''))
predicted.path.c <- lapply(1:ncol(X$Xt.c), function(i) X$Xt.c[,i])
c.rng <- as.numeric(X$c.rng)
# get path data for group
Mice <- path.data[path.data$Condition == grp,-1]
Grp.Mean <- colMeans(Mice)
log.path <- log(Grp.Mean,base=10)
mask <- Grp.Mean != 0 # identify regions with no path --> log = -Inf
fit.by.c <- sapply(predicted.path.c, function(X) cor(log(X,base=10)[mask],log.path[mask]))

p <- ggplot() + geom_col(aes(x=c.rng,y=fit.by.c)) + theme_classic() + xlab('c') + ylab('Fit (r)')
p

ggsave(p,filename=paste(savedir,grp,'FitByCInstability.pdf',sep=''),width=1.5,height=1.5,units='in')
