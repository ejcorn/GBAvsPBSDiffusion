#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)

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
rownames(connectivity.ipsi) <- paste('i',rownames(connectivity.ipsi),sep='') # add i to ipsilateral connections
colnames(connectivity.ipsi) <- paste('i',colnames(connectivity.ipsi),sep='') # add i to ipsilateral connections
rownames(connectivity.contra) <- paste('i',rownames(connectivity.contra),sep='') # add i to ipsilateral connections
colnames(connectivity.contra) <- paste('c',colnames(connectivity.contra),sep='') # add c to contralateral connections

W <- rbind(cbind(connectivity.ipsi,connectivity.contra),cbind(connectivity.contra,connectivity.ipsi))
n.regions <- nrow(W)

# confirm it worked

if(identical(colnames(path.data)[-1],colnames(W))){
  print('Names correctly matched.')
} else{print('Warning!! Names not properly matched.'); break}

# retain indices to reorder like original data variable for plotting on mouse brains
region.names <- colnames(W)
orig.order <- order(match(region.names,names(data)[-1]))

# save data
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))
save(path.data, region.names, orig.order,n.regions, file = paste(savedir,'pathdata.RData',sep=''))

###############################
### Process Snca Expression ###
###############################

Synuclein <- read.csv('data/SncaExpression.csv',header = T,stringsAsFactors = F,check.names = F)

syn.names <- colnames(Synuclein)
syn.names[!syn.names %in% region.names] # names in Snca but not in regions? Yes, delete those
print('regions missing from Snca Expression?')
print(region.names[!region.names %in% syn.names]) # names in regions but not in Snca? No
Synuclein <- Synuclein[,which(syn.names %in% region.names)] # delete regions not in path data
sum(is.na(Synuclein))
Synuclein <- Synuclein[,order(match(colnames(Synuclein),region.names))]
Synuclein <- t(Synuclein)


save(Synuclein,file = paste(savedir,'Snca.RData',sep=''))