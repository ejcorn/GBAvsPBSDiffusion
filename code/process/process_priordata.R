# this script processes the re-annotated data from Henderson et al. 2019 Quantitative Mapping....
# from 1,3,6 MPI

#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/miscfxns.R')

load(file = paste(savedir,'pathdata.RData',sep=''))
Group.Column.new <- path.data[,1]
path.data <- path.data[,-1]	# remove group column, add back in later
prior.data <- read.csv('data/PriorPathData.csv',header = TRUE, check.names = FALSE)
Group.Column.prior <- prior.data[,1] # extracting only MPI, all mice in table are NTG
prior.data <- prior.data[,-(1:2)] # remove time point and group column, add back in later

############################
### Process region names ###
############################

new.path.names <- colnames(path.data) # get names from new path data
prior.path.names <- colnames(prior.data) # get names from prior path data

print('region names in new data not in prior data?')
print(new.path.names[!new.path.names %in% prior.path.names])

print('region names in prior data not in new data?')
print(prior.path.names[!prior.path.names %in% new.path.names])

###############################################
### Exclude regions missing from prior data ###
###############################################

new.missing.mask <- new.path.names %in% prior.path.names	# exclude missing regions
path.data <- path.data[,new.missing.mask]

# get path names again after eliminating regions
new.path.names <- colnames(path.data) # get names from new path data, excluding first 'Group' column

# order prior data to match new data
prior.data <- prior.data[,order(match(prior.path.names,new.path.names))]
unit.test(identical(colnames(prior.data),colnames(path.data)),
	'successfully reordered regions',
	'Regions reordered incorrectly!')

# exclude regions from synuclein
load(file = paste(savedir,'Snca.RData',sep=''))
Synuclein <- matrix(Synuclein[new.missing.mask,],ncol=1)

save(Synuclein,file = paste(savedir,'Snca.RData',sep=''))

# exclude regions from connectome
W <- readMat(paste(savedir,'W.mat',sep=''))
W <- W$W[new.missing.mask,new.missing.mask]
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))

region.names <- colnames(path.data) # get new region names
n.regions <- length(region.names)

# load original GBA data set to preserve original order
orig.data <- read.csv('data/PathData.csv',header = TRUE, check.names = FALSE)
orig.names <- colnames(orig.data)[-1]	# remove Condition column name
orig.order <- order(match(region.names,orig.names))

# confirm orig order maps region.names back to names in orig.names that are also in region.names
unit.test(identical(region.names[orig.order],orig.names[orig.names %in% region.names]),
	'successfully reordered regions',
	'Regions reordered incorrectly!')

# normalize new path data to 1 month from prior data
# normalize separately for CBE or PBS
for(grp in params$grps){
	path.data[Group.Column.new == grp,] <- path.data[Group.Column.new == grp,]*mean(as.matrix(path.data[Group.Column.new == grp,]))/mean(as.matrix(prior.data[Group.Column.prior==1,]))
}


# add back in group/time point columns and save data
path.data <- cbind(Condition=Group.Column.new,path.data)
prior.data <- cbind(MPI=Group.Column.prior,prior.data)

save(path.data, prior.data, region.names, orig.order,orig.names,n.regions, file = paste(savedir,'pathdata.RData',sep=''))
