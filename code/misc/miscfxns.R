unit.test <- function(test,pass,fail){
	# test: logical expression to evaluate
	# pass: error to print if true
	# fail: error to print if false
	if(test){
		print(pass)
	} else{print(fail)}

}

fisher.r.to.z <- function(r){
  r <- 0.5*(log(1+r) - log(1-r))
  return(r)
}