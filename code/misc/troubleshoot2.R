# this script is the same as analyzespread.R
# only it fits time constant using prior data for 3 and 6 months

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params')))
grp = 'DPBS'
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel_priorfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months
tp <- c(1,3,6)
Mice.prior <- lapply(tp[2:3], function(M) prior.data[prior.data$MPI == M,-1])

# join pathology measurements
Mice <- c(list(Mice),Mice.prior)
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

load(paste(params$opdir,'processed/Lout.RData',sep=''))

# Fit time scaling parameter on average of all mice
c.rng <- seq(0.01,1,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
c.Grp <- c.fit.prior(log.path,L.out,Xo,tp,c.rng)

Xt <- predict.Lout(L.out,Xo,c.Grp,t=seq(0,6,length.out = 100))
Xt.name <- do.call('rbind', lapply(region.names, function(R) rep(R,100)))
Xt.t <- do.call('rbind', lapply(region.names, function(R) seq(0,1,length.out = 100)))
df <- data.frame(Xt = as.vector(Xt), grp = as.vector(Xt.name), t = as.vector(Xt.t))
p <- ggplot(data=df) + geom_line(aes(x=t,y=Xt,color=grp)) + theme(legend.position = 'none')
p