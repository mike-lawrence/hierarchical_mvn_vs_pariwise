proper_level_order = function(x){
	(
		x
		%>% unique()
		%>% stringr::str_split(
			string = .
			, pattern = stringr::fixed('.')
			, simplify = T
		)
		%>% {function(z){
			for(col in 1:ncol(z)){
				orig_warn = options(warn=-1)
				tmp = as.numeric(z[,col])
				options(orig_warn)
				if(!any(is.na(tmp))){
					z[,col] = stringr::str_pad(tmp,nchar(max(tmp)),pad='0')
				}
			}
			return(apply(z,1,paste0,collapse='.'))
		}}()
		%>% sort()
		%>% {function(z){
			while(any(stringr::str_detect(z,stringr::fixed('.0')))){
				z %<>% stringr::str_replace(stringr::fixed('.0'),'.')
			}
			return(z)
		}}()
	)
}
