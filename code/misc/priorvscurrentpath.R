# this script compares prior data and current dataset

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

# get mean pathology for only GBA time point, 1 month
Mice <- path.data[path.data$Condition == grp,-1]

# get mean pathology from prior data set 
tp <- c(1,3,6)
Mice.prior <- lapply(tp, function(M) prior.data[prior.data$MPI == M,-1])

# compare prior dataset mean pathology to current dataset mean pathology
prior.means <- lapply(Mice.prior, colMeans)
current.means <- colMeans(Mice)
cor(do.call('cbind',prior.means),current.means)

# compare individual mice between prior and current data set
prior.curr.mice.r <- lapply(1:length(tp), function(M) cor(t(Mice),t(Mice.prior[[M]])))
prior.curr.mice.r <- lapply(prior.curr.mice.r, as.vector)
month <- lapply(1:length(tp), function(M) rep(paste('Month',tp[M]),length(prior.curr.mice.r[[M]])))
df <- data.frame(r=unlist(prior.curr.mice.r),month = unlist(month))

month.means <- sapply(prior.curr.mice.r,mean)
month.sem <- sapply(prior.curr.mice.r,function(x) sd(x)/sqrt(length(x)))
m.names <- sapply(tp, function(m) paste('Month',m))

p <- ggplot() + geom_violin(data=df,aes(x=month,y=r,fill=month)) +
	geom_errorbar(aes(x=m.names,ymin= month.means - 2*month.sem,ymax = month.means + 2*month.sem)) +
	geom_point(aes(x=m.names,y=month.means)) + scale_y_continuous(limits=c(0,1))