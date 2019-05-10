# this script fits time constant using prior data for 3 and 6 months

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'synsweep/',sep='')
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

# loop through different minimum synuclein values
# for scaling W using W = diag(S) %*% W where S is vector of synuclein expression
# min value of 1 means lowest synuclein expression is unchanged
# mean value of 1 means regions with high syn expression increase spread, regions with low decrease
# compute fit and remaining variance explained by synuclein

new.min.rng <- 10^c(-10:-1,1) # specify range of minimum values for synuclein weighting
n.mins <- length(new.min.rng) # number of values to test
# store fit and remaining variance in pathology explained by synuclein
path.fit <- matrix(NA,nrow=n.mins)
syn.var <- matrix(NA,nrow=n.mins)

for(N in 1:n.mins){
	
	new.min.N <- new.min.rng[N]
	print(paste(N,'of',n.mins,'scalings'))
	Syn.norm <- norm.synuclein(S=Synuclein,new.min=new.min.N) # shift and scale synuclein as described above varying new.min
	L.out <- get.Lout(W=W,S=Syn.norm) # recompute out-degree Laplacian
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
	df <- data.frame(path = log.path, Xt = Xt.Grp, Syn = Syn.norm)
	# exclude regions with 0 pathology at each time point for purposes of computing fit
	mask <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
	print(paste(sum(mask),'regions')) # number of regions left after exclusion
	df <- df[mask,] 
	c.test <- cor.test(df$path,df$Xt)
	
	path.fit[N] <- c.test$estimate	# store model fit
	model.resid <- residuals(lm(df$path~df$Xt)) # residuals of Syn-weighted connectivity based prediction
	syn.var[N] <- cor(model.resid,df$Syn) # variance in residuals explained by Synuclein

}

df.plt <- data.frame(x=log(new.min.rng,base=10), Diff.Model = path.fit, Synuclein = syn.var)

# correlation values
p <- ggplot(data=df.plt) + 
	geom_line(aes(x=x,y=Diff.Model),color='light blue') + 
	geom_point(aes(x=x,y=Diff.Model),color='blue') + 
	geom_line(aes(x=x,y=Synuclein),color='darksalmon') +
	geom_point(aes(x=x,y=Synuclein),color='red') +
	annotate(geom='text',x=min(df.plt$x),y=0.9,label='Diffusion Model',color='light blue',hjust=0)+
	annotate(geom='text',x=min(df.plt$x),y=0.7,label='Snca ~ Model Residuals',color='darksalmon',hjust=0)+
	scale_y_continuous(limits=c(0,1)) + ggtitle('') + 
	xlab('log(Min(Synuclein))') + ylab('r')
p

ggsave(p,filename = paste(savedir,grp,'ModelVsSynCorrelationWithPathology.pdf',sep=''),units = 'in',height = 3,width = 4)


# plot variance explained in pathology

df.plt$rsq <- path.fit^2 # variance explained by model
df.plt$syn.rsq <- (1 - df.plt$rsq)*syn.var^2 # amount of remaining variance explained by synuclein

p <- ggplot(data=df.plt) + 
	geom_line(aes(x=x,y=rsq),color='light blue') + 
	geom_point(aes(x=x,y=rsq),color='blue') + 
	geom_line(aes(x=x,y=syn.rsq),color='darksalmon') +
	geom_point(aes(x=x,y=syn.rsq),color='red') +
	annotate(geom='text',x=min(df.plt$x),y=0.7,label='Diffusion Model',color='light blue',hjust=0)+
	annotate(geom='text',x=min(df.plt$x),y=0.5,label='Snca Expression',color='darksalmon',hjust=0)+
	scale_y_continuous(limits=c(0,1)) + ggtitle('Variance in Pathology Explained') + xlab('log(Min(Synuclein))') + ylab('R^2')
p

ggsave(p,filename = paste(savedir,grp,'ModelVsSynR^2Pathology.pdf',sep=''),units = 'in',height = 3,width = 4)
