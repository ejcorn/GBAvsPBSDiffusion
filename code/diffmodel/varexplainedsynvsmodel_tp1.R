# this script fits time constant using only GBA data for 1 month

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel_priorfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W

# get mean pathology for only GBA time point, 1 month
tp <- 1
Mice <- path.data[path.data$Condition == grp,-1]
Grp.mean <- colMeans(Mice,na.rm = T)
log.path <- log(Grp.mean,base=10)
Xo <- get.Xo(region.names,'iCP') # seed pathology in iCPu

# exclude regions with missing synuclein data
# Syn.mask <- !is.na(Synuclein)
# log.path <- log.path[Syn.mask]
# W <- W[Syn.mask,Syn.mask]
# Synuclein <- Synuclein[Syn.mask,]
# Xo <- Xo[Syn.mask]

#S <- norm.synuclein(Synuclein,10^-2)
S <- rep(1,nrow(W))
L.out <- get.Lout(W=W,S=S) # recompute out-degree Laplacian using no synuclein weighting
# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
c.Grp <- c.fit(log.path,L.out,Xo,c.rng)

# for the remainder of the script, only evaluating path at 1 month
# so set log path equal to only 1 month

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
