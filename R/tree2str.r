# tree2str.r #############################################################################################################
# FUNCTION:             		  	DESCRIPTION:
#  get.params						Prints the parameter values.
#  .read.params 					Reads the parameter values from the tree. (Internal function)
#  tree2str							Prints the structure of HAC as string.
#  .allocate.all					Asscosiated with .allocate.one.with.theta and .allocate.one.with.theta. (Internal function)
#  .allocate.one.with.theta			Constructs the string, if theta = TRUE. (Internal function)
#  .allocate.one.without.theta		Constructs the string, if theta = FALSE. (Internal function)
#  ..allocate.all.for.plot          Produces the subsripts for the plot.
##########################################################################################################################

get.params = function(hac, sort.v = FALSE, ...){
		res = numeric(1)
	
    if(hac$type != GAUSS){
    	res = .read.params(hac$tree)
    }else{ 
    	if(hac$type == GAUSS){
        res = "No tree for GAUSS models."
    }}
    if(sort.v == FALSE){res}else{sort(res, ...)}
}

#-------------------------------------------------------------------------------------------------------------------------------

.read.params = function(tree){
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree)
	s = sapply(tree, is.list)
		
		if(any(s)==TRUE){
			res = c(tree[[n]], unlist(sapply(X = tree[s], FUN = .read.params)))
		}else{
			res = tree[[n]]}
	res
}

#-------------------------------------------------------------------------------------------------------------------------------

tree2str = function(hac, theta = TRUE, digits = 2){
	res = character(1)
	
    if(hac$type != GAUSS){
        res = .allocate.all(hac$tree, theta, digits)
    }else{ 
    	if(hac$type == GAUSS){
        res = "No tree for GAUSS models"
    }}
    res
}

#-------------------------------------------------------------------------------------------------------------------------------

.allocate.all = function(tree, theta, digits){
n = length(tree); x = character(1)
	
		if(theta == TRUE){
			for(i in 1:n){
				if(i == 1){
					x = paste("(", .allocate.one.with.theta(tree[[i]], digits = digits, theta = theta), sep = "")
			}else{
				if((i > 1) & (i < n)){
					x = paste(x, ".", .allocate.one.with.theta(tree[[i]], digits = digits, theta = theta), sep = "")
			}else{
				if(i == n){
					x = paste(x, ")_{", .allocate.one.with.theta(tree[[i]], digits = digits, theta = theta),"}", sep = "")
			}}}}
	}else{
		if(theta == FALSE){
			for(i in 1:(n-1)){
				if(i == 1){
					x = paste("(", .allocate.one.without.theta(tree[[i]], theta = theta), sep = "")
			}else{
				if((i > 1) & (i < n)){
					x = paste(x, ".", .allocate.one.without.theta(tree[[i]], theta = theta), sep = "")
			}}}}
					x = paste(x, ")", sep = "")
			}
	return(x)
}

#-------------------------------------------------------------------------------------------------------------------------------

.allocate.one.with.theta = function(element, digits, theta){
		if(class(element)=="character"){
			d = length(element)
			if(d ==1){
				element
		}else{
			if(d > 1){
				paste(element, sep = "", collapse = ".")
		}}
	}else{
		if(class(element)=="list"){
			.allocate.all(element, theta = theta, digits = digits)
	}else{
		if(class(element)=="numeric"){
			round(element, digits = digits)
		}}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.allocate.one.without.theta = function(element, theta){
		if(class(element)=="character"){
			d = length(element)
			if(d ==1){
				element
		}else{
			if(d > 1){
				paste(element, sep = "", collapse = ".")
		}}
	}else{
		if(class(element)=="list"){
			return(.allocate.all(element, theta = theta))
	}}
}
