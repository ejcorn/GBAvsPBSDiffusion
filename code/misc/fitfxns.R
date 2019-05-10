get.Xo <- function(region.names,ROI){
  # generate initial "pathology seed" vector, which is all 0's except
  # a 1 in one region determined by character vector ROI
  # ROI should correspond to region.names
  n.regions <- length(region.names)
  Xo <- matrix(0, nrow=n.regions,ncol=1)
  Xo[which(region.names == ROI)] <- 1
  return(Xo)
}


norm.synuclein <- function(S,new.min){
  # S: vector of synuclein expression values
  # new.min: value to specify for new minimum value in synuclein vector
  # this function will set the minimum equal to new.min and the mean equal to 1

  S <- S - min(S,na.rm=T) + new.min # can't disconnect graph
  S <- S/mean(S) # set mean to 1 
  return(S)

}

get.Lout <- function(W,S){
  # W: NxN adjacency matrix, may be asymmetric (if so, connections are from row i to column j)  
  # S: vector of weights to apply to each row of W (here, synuclein expression)
  # The connections directionallity here is opposite the convention for matrix multiplication
  # so we are capturing retrograde spread along connections
  # compute out-degree laplacian for retrograde spread

  n.regions <- nrow(W)
  W <- W * !diag(n.regions) # get rid of diagonal
  W <- diag(as.numeric(S)) %*% W
  #W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

  # Where i is row element and j is column element
  # Wij is a connection from region i to region j
  # convention is the opposite, so without transposing W
  # I am capturing "retrograde" connectivity

  in.deg <- colSums(W)
  out.deg <- rowSums(W)
  L.out <- diag(x = out.deg) - W # outdegree laplacian
  return(L.out)
}

predict.Lout <- function(L,Xo,c,t=1){
  # Generate prediction for given time points using
  # L: Laplacian of adjacency matrix
  # Xo: initial conditions
  # c: time constant
  # t: time points

  Xt <- do.call('cbind', lapply(t, function(t.i) expm(-L*c*t.i)%*%Xo))
  return(Xt)
}

c.fit <- function(log.path,L.out,Xo,c.rng){
  n.regions <- length(log.path)
  # exclusion mask... don't count regions with 0 path
  mask <- log.path != -Inf
  # compute fit at each time point for range of time
  Xt.sweep <- sapply(c.rng, function(c) # no time here, because in this project there's only 1 time pointo 
    cor(log.path[mask],log(predict.Lout(L.out,Xo,c),base=10)[mask]))
  c <- c.rng[which.max(Xt.sweep)] # select c giving max correlation
  print(c)
  return(c)
}

c.fit.r <- function(log.path,L.out,Xo,c.rng){
  n.regions <- length(log.path)
  # exclusion mask... don't count regions with 0 path
  mask <- log.path != -Inf
  # compute fit at each time point for range of time
  Xt.sweep <- sapply(c.rng, function(c) # no time here, because in this project there's only 1 time pointo 
    cor(log.path[mask],log(predict.Lout(L.out,Xo,c),base=10)[mask],use='pairwise.complete.obs'))
  c <- c.rng[which.max(Xt.sweep)] # select c giving max correlation
  r <- max(Xt.sweep)
  print(c)
  return(list(c=c,r=r,Xt.sweep=Xt.sweep)) # return time constant and correlation coefficient
}

c.fit.prior <- function(log.path,L.out,Xo,tp,c.rng){
  # log.path: list of empirically measured path values for each node
  # L.out: Laplacian matrix based on structural connectivity
  # tp: vector of real values for time
  # Xo: vector of initial state of pathology
  # c.rng: values of time constant to test

  n.regions <- length(log.path[[1]])
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  # compute fit at each time point for range of time
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(predict.Lout(L.out,Xo,c,t=t),base=10)[mask[[t]]])))
  Xt.sweep <- Reduce('+',Xt.sweep) / length(tp) # mean fit for each tp
  c <- c.rng[which.max(Xt.sweep)]
  print(c)
  return(c)
}

c.fit.prior.r <- function(log.path,L.out,Xo,tp,c.rng){
  # log.path: list of empirically measured path values for each node
  # L.out: Laplacian matrix based on structural connectivity
  # tp: vector of real values for time
  # Xo: vector of initial state of pathology
  # c.rng: values of time constant to test

  n.regions <- length(log.path[[1]])
  # exclusion mask... don't count regions with 0 path
  mask <- lapply(log.path, function(x) x != -Inf)
  # compute fit at each time point for range of time
  Xt.sweep <- lapply(1:length(tp), function(t) # continuous time
    sapply(c.rng, function(c)
      cor(log.path[[t]][mask[[t]]],log(predict.Lout(L.out,Xo,c,t=t),base=10)[mask[[t]]])))
  r <- sapply(Xt.sweep, function(x) max(x)) # best fit for each time point
  Xt.sweep <- Reduce('+',Xt.sweep) / length(tp) # mean fit for each tp
  c <- c.rng[which.max(Xt.sweep)]
  print(c)
  return(list(c=c,r=r))
}