# constructor.r ##########################################################################################################
# FUNCTION:         DESCRIPTION:
#  .check.par			  Checks, whether the parameter are correctly specified, i.e., ascending ordered from lowest to the highest hierarchical level. (Internal function)
#  .check				    Checks, whether the parameters of two subsequent nodes are correctly specified. (Internal function)
#  hac.full				  Constructs 'hac' objects for fully nested Archimedean Copulae.
#  hac					    Constructs 'hac' objects for arbitrary nested Archimedean Copulae.
#  print.hac     		Prints 'hac' objects.
##########################################################################################################################

.check.par = function(x){
	if((x$type == AC_GUMBEL) || (x$type == HAC_GUMBEL)){
		ober.theta = 1}
	else
	if((x$type = AC_CLAYTON) || (x$type = HAC_CLAYTON)){
		ober.theta = 0}

	L = x$tree
	.check(L = L, theta = ober.theta)
}

#------------------------------------------------------------------------------------------------------------------------

.check = function(L, theta){
	n = length(L)
	if(L[[n]] < theta){
		return(warning("The dependency parameter of the nested AC should be higher than the parameter at the initial node."))
	}else{
	for(i in 1:(n-1)){
		if(class(L[[i]]) == "list"){
			.check(L = L[[i]], theta = L[[n]])
	}}}
}

#------------------------------------------------------------------------------------------------------------------------

hac.full = function(type = HAC_GUMBEL, y, theta){
	n = length(y)
	if(n != (length(theta) + 1)){return(warning("The input arguments does not fit to a fully nested HAC"))}
	
	tree = list(y[n], y[n-1], theta[n-1])
	for(i in (n-2):1){
		tree = list(tree, y[i], theta[i])
	}
	hac(type = type, tree = tree)
}

#--------------------------------------------------------------------------------------------
 
hac = function(type = HAC_GUMBEL, tree = NULL, corr = NULL){
 	object = list(type = type, tree = tree)
 	class(object) = "hac"
    if(type != GAUSS){.check.par(object)}else{object$tree = corr} 
 	object
}

#------------------------------------------------------------------------------------------------------------------------
	
print.hac = function(x, digits = 2, ...){
        cat("Class: hac", "\n", ...)
	if((x$type == 0) | (x$type == 1)){
 		   cat("Generator: Gumbel", "\n", ...)
 		   cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 		}
	if((x$type == 3) | (x$type == 4)){
 		   cat("Generator: Clayton", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 		}
    if(x$type == 5){
        cat("Family: Gaussian", "\n", ...)
        cat(round(x$tree[lower.tri(x$tree)], digits = digits, ...), "\n", ...)
 	}
 }
