plot.Xt <- function(Xt,t){
	# Xt: NxT matrix of nodal value (N) over time (T)
	t <- do.call('rbind',lapply(1:length(Xo), function(i) matrix(t,nrow=1)))
	ROI <- do.call('cbind',lapply(t, function(i) 1:length(Xo)))
	df <- data.frame(y=as.vector(Xt),x=as.vector(t),as.vector(ROI))
	p <- ggplot(df) + geom_line(aes(x=x,y=y,color=ROI))
	p

}