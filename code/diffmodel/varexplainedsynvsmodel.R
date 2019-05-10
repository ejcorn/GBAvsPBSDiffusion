# this script fits time constant using prior data for 3 and 6 months

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel_priorfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months
tp <- c(1,3,6)
Mice.prior <- lapply(tp[2:3], function(M) prior.data[prior.data$MPI == M,-1])

# join pathology measurements
Mice <- c(list(Mice),Mice.prior)
Grp.mean <- lapply(Mice, function(x) colMeans(x,na.rm = T))

S <- norm.synuclein(Synuclein,10^-2)
#S <- rep(1,nrow(W))
L.out <- get.Lout(W=W,S=S) # recompute out-degree Laplacian using no synuclein weighting
# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(x) log(x,base=10))
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu
c.Grp <- c.fit.prior(log.path,L.out,Xo,tp,c.rng)

# for the remainder of the script, only evaluating path at 1 month
# so set log path equal to only 1 month
log.path <- unlist(log.path[1])
tp <- 1

Xt.Grp <- log(predict.Lout(L.out,Xo,c.Grp), base = 10) # predict pathology using connectivity, time constant, and seed
#Xt.Grp <- scale(Xt.Grp,center=T) # *** added 2/15 to visualize whether correlation was real or not
df <- data.frame(path = log.path, Xt = Xt.Grp, Syn = Synuclein)
# exclude regions with 0 pathology at each time point for purposes of computing fit
mask <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
print(paste(sum(mask),'regions')) # number of regions left after exclusion
df <- df[mask,] 
c.test <- cor.test(df$path,df$Xt)
print(paste('diffusion R^2 =',c.test$estimate^2))
c.test.syn <- cor.test(df$path,df$Syn)
print(paste('synuclein R^2 =',c.test.syn$estimate^2))

full.mdl <- lm(path~Xt + Syn,data=df)
print(paste('max possible R^2 =',summary(full.mdl)$r.squared))