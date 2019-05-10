rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion/'
setwd(basedir)
params <- list(basedir=basedir,
               matlab.path='/Applications/MATLAB/R2018b/bin/matlab',
               grps = c('DPBS','CBE'),
               c.min = 0.01, # empirically conservative minimum for time constant
               c.max = 100, # empirically conservative maximum for time constant
               c.n = 100) # number of time constants to try
params$opdir <- paste('GBAvsPBSasyndiffusionCMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#####################
### Load packages ###
#####################

source('code/misc/packages.R')

####################
### Process data ###
####################

source('code/process/process.R')
source('code/process/getLout.R')

################################
### Run base diffusion model ###
################################

for(grp in params$grps){
  source('code/diffmodel/analyzespread.R')
  source('code/diffmodel/seedspec.R')
  source('code/diffmodel/examineseedspec.R')
}

###############################################
### Compare distributions of time constants ###
###############################################

for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEtimeconstbymouse.R')}
source('code/PBSvsCBE/plotPBSvsCBEtimeconstbymouse.R')