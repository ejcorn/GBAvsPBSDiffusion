plot.Xt <- function(Xt,t){
	# Xt: NxT matrix of nodal value (N) over time (T)
	t <- do.call('rbind',lapply(1:length(Xo), function(i) matrix(t,nrow=1)))
	ROI <- do.call('cbind',lapply(t, function(i) 1:length(Xo)))
	df <- data.frame(y=as.vector(Xt),x=as.vector(t),as.vector(ROI))
	p <- ggplot(df) + geom_line(aes(x=x,y=y,color=ROI))
	p

}

library(reshape2)
imagesc <- function(m){
	# plot matrix m as a heatmap for quick visualization
	df <- melt(as.matrix(m))
	p <- ggplot(df) + geom_tile(aes(x=Var1,y=Var2, fill=value)) + 
		scale_y_discrete(limits=rev(unique(df$Var1))) + # flip y-axis so region 1 is at the top of the plot, 
		theme(axis.text.x = element_text(angle=90))
	return(p)
}