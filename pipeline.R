rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion/'
setwd(basedir)
params <- list(basedir=basedir,               
               grps = c('DPBS','CBE'),
               c.min = 0.01, # empirically conservative minimum for time constant
               c.max = 10, # empirically conservative maximum for time constant
               c.n = 100) # number of time constants to try
source('code/misc/miscfxns.R')
params$source.save <- source.save
params$opdir <- paste('GBAvsPBSasyndiffusion112519CMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#################################################
### Load packages & create output directories ###
#################################################

source('code/misc/packages.R')
source('code/misc/directories.R')

####################
### Process data ###
####################

source('code/process/process.R')
source('code/process/getLout.R')

################################
### Run base diffusion model ###
################################

for(grp in params$grps){
  # save all stats output
  params$source.save('code/diffmodel/analyzespread.R',paste(params$opdir,'diffmodel/',grp,'spread.log',sep=''))
  source('code/diffmodel/seedspec.R')
  source('code/diffmodel/plotseedspec.R')
  # save all stats output
  params$source.save('code/diffmodel/examineseedspec.R',paste(params$opdir,'diffmodel/seedspec/',grp,'connsimaltseed.log',sep=''))
}



###############################################
### Compare distributions of time constants ###
###############################################

source('code/PBSvsCBE/predict_crng.R')
for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEtimeconstbymouse.R')}
source('code/PBSvsCBE/plotPBSvsCBEtimeconstbymouse.R')
for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEfitbycdistribution.R')}
source('code/PBSvsCBE/plotPBSvsCBEfitbycdistribution.R')

for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEtimeconst_boot.R')}
source('code/PBSvsCBE/plotPBSvsCBEtimeconst.R')
