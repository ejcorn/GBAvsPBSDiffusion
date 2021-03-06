#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/miscfxns.R')

data <- read.csv('data/PathData.csv',header = TRUE, check.names = FALSE)
connectivity.ipsi <- read.csv('data/Connectome_Ipsi.csv',row.names = 1,header = TRUE,check.names = FALSE)
connectivity.contra <- read.csv('data/Connectome_Contra.csv',row.names = 1,header = TRUE, check.names = FALSE)

############################
### Process region names ###
############################

conn.names.ipsi <- colnames(connectivity.ipsi)
conn.names.contra <- colnames(connectivity.contra)

# checks 
if(identical(colnames(connectivity.contra),rownames(connectivity.contra))){
  print('contra connectivity colnames and rownames equal')}
if(identical(colnames(connectivity.ipsi),rownames(connectivity.ipsi))){
  print('ipsi connectivity colnames and rownames equal')}

path.names <- colnames(data)[-1] # get names from path data, excluding first 'Group' column
#find ipsi/contralateral regions as denoted by the c/i as the first character
cInd <- which(sapply(path.names, function(X) substr(X,1,1) == 'c'))  # index of regions contra to injection site
iInd <- which(sapply(path.names, function(X) substr(X,1,1) == 'i'))  # index of regions ipsi to injection site
# extract ipsi/contra regions
path.names.contra <- unname(sapply(cInd, function(X) substr(path.names[X],2,nchar(path.names[X]))))
path.names.ipsi <- unname(sapply(iInd, function(X) substr(path.names[X],2,nchar(path.names[X]))))

# do regions match between connectome and path data?
if(identical(path.names.ipsi,path.names.contra)){
  print('Ipsi and contra path names identical')} # should be TRUE
if(identical(conn.names.contra,conn.names.ipsi)){
  print('Ipsi and contra connectivity names identical')} # should be TRUE

print('extra names in path?')
print(path.names.ipsi[!path.names.ipsi %in% conn.names.ipsi]) # any extra names in path.names?
print('extra names in connectivity?')
print(conn.names.ipsi[!conn.names.ipsi %in% path.names.ipsi]) # any extra names in conn.names?

# order path and connectivity by same names

path.data <- data[,-1] # get path data without group column
path.data.ipsi <- path.data[,iInd]
path.data.contra <- path.data[,cInd]
# concatenate ipsi then contra in an order that matches connectome
# and add back in condition column
path.data <- cbind(data[,1],path.data.ipsi[,order(match(path.names.ipsi,conn.names.ipsi))],
                   path.data.contra[,order(match(path.names.contra,conn.names.contra))]) 

colnames(path.data)[1] <- 'Condition'

# tile matrix such that sources are rows, columns are targets (see Oh et al. 2014 Fig 4)

W <- rbind(cbind(connectivity.ipsi,connectivity.contra),cbind(connectivity.contra,connectivity.ipsi))
rownames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral connections
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral connections
colnames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral connections
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral connections

n.regions <- nrow(W)

# confirm it worked

if(identical(colnames(path.data)[-1],colnames(W))){
  print('Names correctly matched.')
} else{print('Warning!! Names not properly matched.'); break}

# retain indices to reorder like original data variable for plotting on mouse brains
region.names <- colnames(W)
orig.names <- names(data)[-1]
orig.order <- order(match(region.names,orig.names))

unit.test(identical(region.names[orig.order],orig.names),'orig.order correct','orig.order INCORRECT')

# save data
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))
save(path.data, region.names, orig.order,orig.names, n.regions, file = paste(savedir,'pathdata.RData',sep=''))

###############################
### Process Snca Expression ###
###############################

Synuclein <- read.csv('data/SncaExpression.csv',header = T,stringsAsFactors = F,check.names = F)

syn.names <- colnames(Synuclein)
syn.names[!syn.names %in% region.names] # names in Snca but not in regions? Yes, delete those
print('regions missing from Snca Expression?')
print(region.names[!region.names %in% syn.names]) # names in regions but not in Snca? No
Synuclein <- Synuclein[,which(syn.names %in% region.names)] # delete regions not in path data
unit.test(identical(orig.names,colnames(Synuclein)),'names check out','SOMETHING IS WRONG')
sum(is.na(Synuclein))
orig.names.syn <- orig.names[!is.na(Synuclein)] # get names in original order, minus those missing Snca data
Synuclein <- Synuclein[,order(match(colnames(Synuclein),region.names))]
Synuclein <- t(Synuclein)
unit.test(identical(rownames(Synuclein),region.names),'synuclein names check out','SOMETHING IS WRONG')

# retain indices to reorder like original data variable for plotting on mouse brains, only using regions with available Snca expression data
region.names.syn <- region.names[!is.na(Synuclein)]
orig.order.syn <- order(match(region.names.syn,orig.names.syn))
# check if we correctly reorder connectome ordered names, minus regions with no syn data, match up with original order from input data, minus regions with no syn data
unit.test(identical(region.names.syn[orig.order.syn],orig.names.syn),'orig.order.syn correct','orig.order.syn INCORRECT')

save(Synuclein,orig.names.syn,orig.order.syn,file = paste(savedir,'Snca.RData',sep=''))
