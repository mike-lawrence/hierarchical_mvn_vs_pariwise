guess_if_intercept = function(x){
	x = as.data.frame(x)
	is_intercept = rep(F,ncol(x))
	for(i in 1:ncol(x)){
		vals = sort(unique(x[,i]))
		if(length(vals)==1){
			if(vals==1){
				is_intercept[i] = T
			}
		}else{
			if(length(vals)==2){
				if(all(vals==c(0,1))){
					is_intercept[i] = T
				}
			}
		}
	}
	return(is_intercept)
}
