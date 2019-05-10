#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression

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

# for the remainder of the script, only evaluating path at 1 month
# so set log path equal to only 1 month

Xt.Grp <- log(predict.Lout(L.out,Xo,c.Grp), base = 10) # predict pathology using connectivity, time constant, and seed
Xt.Grp <- scale(Xt.Grp,center=T)
df <- data.frame(path = log.path, Xt = Xt.Grp, Syn = as.numeric(Synuclein))
# exclude regions with 0 pathology at each time point for purposes of computing fit
mask <- df$path != -Inf
print(paste(sum(mask),'regions')) # number of regions left after exclusion
df <- df[mask,] 
c.test <- cor.test(df$path,df$Xt)
print(paste(grp,'p =',c.test$p.value))
p <- ggplot(df,aes(x=Xt,y=path)) + geom_smooth(color = '#007257',method ='lm',size=1) + geom_point(color = '#007257',size = 1,alpha=0.6,stroke=0) +
    annotate(geom='text',x=max(df$Xt) - 0.2*diff(range(df$Xt)),y=min(df$path) + 0.1,label = paste('r =',signif(c.test$estimate,2)),size=2.5) +
    #scale_x_continuous(limits=c(min(df$Xt),max(df$Xt))) +
    theme_classic() + xlab('log(Predicted)') + ylab('log(Path.)') + ggtitle(grp) +
    theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8))

ggsave(p,filename = paste(savedir,grp,'predictedpathcontinuous.pdf',sep=''),units = 'in',height = 1.5,width = 1.5)

# define vulnerability as pathology with connectivity and synuclein regressed out 
vulnerability <- residuals(lm(path~Xt+Syn,data=df))
# *Eli* indicates that it regions are in the order of Oh et al. connectome, which is used throughout code
write.csv(as.data.frame(vulnerability),paste(savedir,grp,'vulnerabilityEli.csv',sep=''),row.names = F)
write.csv(as.data.frame(Xt.Grp),paste(savedir,grp,'predictedpathEli.csv',sep=''),row.names = F)
v.names <- region.names

# reorder regions for plotting on mouse brains (relative to original data)
vulnerability.df <- data.frame(vuln=vulnerability[orig.order],pred=Xt.Grp[orig.order])
rownames(vulnerability.df) <- region.names[orig.order]
write.csv(vulnerability.df,paste(savedir,grp,'vulnerabilitypredicted.csv',sep=''),row.names=FALSE)
