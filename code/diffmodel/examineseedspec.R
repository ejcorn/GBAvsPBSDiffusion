#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/seedspec/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set for 3 and 6 months

Grp.mean <- colMeans(Mice,na.rm = T)
# Fit time scaling parameter on average of all mice
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W

######################################################
### Examine seed regions with better fits than CPu ###
######################################################

load(file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
W <- W * !diag(n.regions) # get rid of diagonal

# extract fits using seeds that aren't the iCP or cCP
# excluding cCP because projection similarity is identical
alt.seeds <- which(rownames(seed.fits) != 'iCP' & rownames(seed.fits) != 'cCP')
iCP <- which(rownames(seed.fits) == 'iCP')
null.cors <- seed.fits[alt.seeds,]

# compute similarity of connectivity to iCPU

conn.sim.out <- sapply(alt.seeds, function(R) cor(W[R,],W[iCP,]))
conn.sim.in <- sapply(alt.seeds, function(R) cor(W[,R],W[,iCP]))
#conn.sim.out <- fisher.r.to.z(conn.sim.out)
#conn.sim.in <- fisher.r.to.z(conn.sim.in)

# plot

df <- data.frame(conn = conn.sim.in,nullcor = null.cors)
m <- gam(nullcor ~ s(conn,k=3),data = df,method = 'REML')
p.val <- signif(summary(m)$s.table[,'p-value'],2)
print(grp)
print('conn.sim.in')
print(summary(m))
print(paste('p =',p.val))
r.sq <- signif(summary(m)$r.sq,2)
lb <- paste('R^2 ==',r.sq)
print(lb)
# "REML and ML are less prone to local minima and may be preferred" 
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html
p.spec <- ggplot(data=df,aes(x=conn,y=nullcor)) + geom_point(color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1) + 
  geom_smooth(method='gam',formula=y~s(x,k=3),color ='#5F4B8B') +
  annotate(geom='text',label=lb,x=max(df$conn),y=min(df$nullcor),hjust=1,vjust=0,parse=T,size=2) +
  xlab('In-Projection Similarity \n to iCPu (r)') + ylab('Fit to iCPu Injection (r)') + ggtitle(grp) +
  theme_classic() + theme(text=element_text(size=6),plot.title = element_text(size=8,hjust=0.5)) #+
  #theme(axis.title.x = element_text(size=8))

ggsave(p.spec,filename = paste(savedir,grp,'AlternateSeedFitByInProjectionSimilarity.pdf',sep=''),units = 'in',height = 1.5,width = 1.5)

df <- data.frame(conn = conn.sim.out,nullcor = null.cors)
m <- gam(nullcor ~ s(conn,k=3),data = df,method = 'REML')
p.val <- signif(summary(m)$s.table[,'p-value'],2)
print(grp)
print('conn.sim.out')
print(summary(m))
print(paste('p =',p.val))
r.sq <- signif(summary(m)$r.sq,2)
lb <- paste('R^2 ==',r.sq)
print(lb)
# "REML and ML are less prone to local minima and may be preferred" 
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html
p.spec <- ggplot(data=df,aes(x=conn,y=nullcor)) + geom_point(color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1) + 
  geom_smooth(method='gam',formula=y~s(x,k=3),color ='#5F4B8B') +
  annotate(geom='text',label=lb,x=max(df$conn),y=min(df$nullcor),hjust=1,vjust=0,parse=T,size=2) +
  xlab('Out-Projection Similarity \n to iCPu (r)') + ylab('Fit to iCPu Injection (r)') + ggtitle(grp) +
  theme_classic() + theme(text=element_text(size=6),plot.title = element_text(size=8,hjust=0.5)) #+
  #theme(axis.title.x = element_text(size=8))

ggsave(p.spec,filename = paste(savedir,grp,'AlternateSeedFitByOutProjectionSimilarity.pdf',sep=''),units = 'in',height = 1.5,width = 1.5)
