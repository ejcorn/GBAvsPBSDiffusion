#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/',sep='')
dir.create(savedir,recursive=T)
#hdir.create(paste(savedir,'roilevel',sep=''),recursive = T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for only time point, 1 month
tp <- 1
Mice <- path.data[path.data$Condition == grp,-1]
Grp.mean <- colMeans(Mice,na.rm = T)

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- log(Grp.mean,base=10)
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
c.Grp <- c.fit(log.path,L.out,Xo,c.rng)

# plot model time course
t.max <- 25
Xt <- predict.Lout(L.out,Xo,1,t=seq(0,t.max,length.out = 100))
Xt.name <- do.call('rbind', lapply(region.names, function(R) rep(R,100)))
Xt.t <- do.call('rbind', lapply(region.names, function(R) seq(0,t.max,length.out = 100)))
df <- data.frame(Xt = as.vector(Xt), grp = as.vector(Xt.name), t = as.vector(Xt.t))
p <- ggplot(data=df) + geom_line(aes(x=t,y=Xt,color=grp),size=0.25) + ylab('x(t)') +
	theme(legend.position = 'none') + theme(text=element_text(size=8))

#ggsave(p,filename = paste(savedir,grp,'TroubleShoot.pdf',sep=''),units = 'in',height = 2,width = 2)
ggsave(p,filename = paste(savedir,grp,'TroubleShoot.png',sep=''),units = 'in',height = 2,width = 2)
