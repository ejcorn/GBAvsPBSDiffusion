rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'PBSvsCBE/',sep='')

tf = 1 # this variable is an irrelevant placehoder

load(paste(savedir,'DPBSSyntimeconstantsbymouse_priorfitTF',tf,'.RData',sep=''))
c.train.PBS <- c.train.Grp
r.PBS <- r.Grp
load(paste(savedir,'CBESyntimeconstantsbymouse_priorfitTF',tf,'.RData',sep=''))
c.train.CBE <- c.train.Grp
r.CBE <- r.Grp
n.reps <- length(c.train.Grp)

df <- data.frame(y=c(c.train.PBS,c.train.CBE),
                 x=c(rep('PBS',n.reps),rep('CBE',n.reps)))

w.t <- wilcox.test(x=c.train.PBS,y=c.train.CBE,conf.int=T)


p <- ggplot(df,aes(y=y,x=x)) + geom_jitter(alpha=0.5,size=1,stroke=0) + ylab('Time Constant') + theme_classic()+
  theme(text=element_text(size=8)) + xlab('') + geom_text(aes(x='CBE',y= 0.35,label = paste('p =',signif(w.t$p.value,2))),size=2)
p
ggsave(plot = p,filename = paste(savedir,'SynTCjitterbymouse_priorfitTF',tf,'.pdf',sep=''), height=2,width=2,units='in')

low.c.fits <- r.PBS[which(c.train.PBS < 0.4)]
hi.c.fits <- r.PBS[which(c.train.PBS > 0.4)]

df = data.frame(PBS=c.train.PBS,CBE=c.train.CBE)
write.csv(x=df,file=paste(savedir,'SynTCbymouse_priorfit.csv',sep=''))