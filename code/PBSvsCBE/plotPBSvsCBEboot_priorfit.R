rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'PBSvsCBE/',sep='')

tf = 'SRS'

load(paste(savedir,'DPBSSyntimeconstants_priorfitTF',tf,'.RData',sep=''))
c.train.PBS <- c.train.Grp
r.PBS <- r.Grp
load(paste(savedir,'CBESyntimeconstants_priorfitTF',tf,'.RData',sep=''))
c.train.CBE <- c.train.Grp
r.CBE <- r.Grp
n.reps <- length(c.train.Grp)

# treat each bootstrapped sample as independent
# between all bootstrapped samples, what percentage is PBS n.eq. CBE
# i.e. a two-tailed test for difference in means
boot.PBS <- matrix(c.train.PBS,nrow=n.reps,ncol=n.reps)
boot.CBE <- t(matrix(c.train.CBE,nrow=n.reps,ncol=n.reps))
boot.diff.means <- as.vector(boot.PBS - boot.CBE)
boot.diff.means <- c.train.PBS - c.train.CBE
# 95% CI
CI95 <- quantile(boot.diff.means,c(0.025,0.975))

df <- data.frame(y=boot.diff.means)    


p <- ggplot(df) + geom_histogram(aes(x=y),alpha=0.5,size=1) + ylab('Count') + 
	annotate(geom='text',x=-0.2,y= 15,label = paste('95% CI: ',signif(CI95[1],2),'-',signif(CI95[2],2),sep=''),size=2,hjust=0) +
	theme_classic() + theme(text=element_text(size=8)) + xlab('Time Constant, PBS - CBE') 
p
ggsave(plot = p,filename = paste(savedir,'SynTCbootdiffjitter_priorfitTF',tf,'.pdf',sep=''), height=2,width=2,units='in')

