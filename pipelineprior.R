rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/GBAvsPBSDiffusion/'
setwd(basedir)
params <- list(basedir=basedir,
               matlab.path='/Applications/MATLAB/R2018b/bin/matlab',
               grps = c('DPBS','CBE'),
               c.min = 0.01, # empirically conservative minimum for time constant
               c.max = 1, # empirically conservative maximum for time constant
               c.n = 100) # empirically conservative number of time constants to try
params$opdir <- paste('GBAvsPBSasyndiffusionpriorfitCMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#####################
### Load packages ###
#####################

source('code/misc/packages.R')

####################
### Process data ###
####################

source('code/process/process.R')
source('code/process/process_priordata.R')
source('code/process/getLout.R')

###############################################
### Troubleshoot normalization of synuclein ###
###############################################
grp <- params$grps[1]
source('code/diffmodel/sweepsynnorm.R')

################################
### Run base diffusion model ###
################################

#for(grp in params$grps){source('code/diffmodel/analyzespread.R')}
#for(grp in params$grps){source('code/diffmodel/analyzespread_prior.R')}
grp <- params$grps[1]
  source('code/diffmodel/analyzespread_prior.R')
  source('code/diffmodel/seedspec_prior.R')

###############################################
### Compare distributions of time constants ###
###############################################

#for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEsyntimeconst.R')}
#source('code/PBSvsCBE/plotPBSvsCBEtimeconst.R')

# for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEsyntimeconst_priorfit.R')}
# for(grp in params$grp){source('code/PBSvsCBE/PBSvsCBEsyntimeconstbymouse_priorfit.R')}
# source('code/PBSvsCBE/plotPBSvsCBEtimeconst_priorfit.R')
