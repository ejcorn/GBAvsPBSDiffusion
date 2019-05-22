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

load(file=paste(savedir,'CBEFitsByCByMouse.RData',sep=''))
df.CBE <- df[seq(1,nrow(df),length.out=100),]
load(file=paste(savedir,'DPBSFitsByCByMouse.RData',sep=''))
df.DPBS <- df[seq(1,nrow(df),length.out=100),]
df <- rbind(df.CBE,df.DPBS)

p <- ggplot(df) + geom_errorbar(aes(x=x,ymin=m-s,ymax=m+s,color=grp),size=0.1) + theme_classic() +
	xlab('Time Constant (c)') + ylab('Fit (r)') + theme(text=element_text(size=8))
p

ggsave(plot = p,filename = paste(savedir,'CBEvsDPBSFitsByTimeConstant.pdf',sep=''), height=2,width=4,units='in')