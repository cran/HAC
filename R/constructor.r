# constructor.r ##########################################################################################################
# FUNCTION:         DESCRIPTION:
#  .check.par		Checks, whether the parameter are correctly specified, i.e., ascending ordered from lowest to the highest hierarchical level. (Internal function)
#  .check			Checks, whether the parameters of two subsequent nodes are correctly specified. (Internal function)
#  hac.full			Constructs 'hac' objects for fully nested Archimedean Copulae.
#  hac				Constructs 'hac' objects for arbitrary nested Archimedean Copulae.
#  print.hac     	Prints 'hac' objects.
#  hac2nacopula     Converts an 'hac' object into a 'nacopula' object.
#  nacopula2hac     Converts a 'nacopula' object into an 'hac' object.
##########################################################################################################################

.check.par = function(x){
	if((x$type == AC_GUMBEL) || (x$type == HAC_GUMBEL) || (x$type == AC_JOE) || (x$type == HAC_JOE)){
		ober.theta = 1}
	else
	if((x$type = AC_CLAYTON) || (x$type = HAC_CLAYTON) || (x$type = AC_FRANK) || (x$type = HAC_FRANK)){
		ober.theta = 1e-10}
	if((x$type = AC_AMH) || (x$type = HAC_AMH)){
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
 
hac = function(type = HAC_GUMBEL, tree = NULL){
 	object = list(type = type, tree = tree)
 	class(object) = "hac"
 	if(((type == AC_AMH) || (type == HAC_AMH)) && (max(.read.params(tree)) >= 1)){return(warning("The largest parameter of the Ali-Mikhail-Haq family must be < 1"))}
    .check.par(object);
 	object
}

#------------------------------------------------------------------------------------------------------------------------
	
print.hac = function(x, digits = 2, ...){
        cat("Class: hac", "\n", ...)
	if((x$type == 1) | (x$type == 2)){
 		   cat("Generator: Gumbel", "\n", ...)
 		   cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
	if((x$type == 3) | (x$type == 4)){
 		   cat("Generator: Clayton", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 5) | (x$type == 6)){
 		   cat("Generator: Frank", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 7) | (x$type == 8)){
 		   cat("Generator: Joe", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 	if((x$type == 9) | (x$type == 10)){
 		   cat("Generator: Ali-Mikhail-Haq", "\n", ...)
	     cat(tree2str(x, theta = TRUE, digits = digits),  "\n", ...)
 	}
 }

#------------------------------------------------------------------------------------------------------------------------

hac2nacopula = function(x){
	.family = character(1)
	.names = .get.leaves(x$tree)
	if((x$type == HAC_GUMBEL) | (x$type == AC_GUMBEL)){
		.family = "G"
	}else
	if((x$type == HAC_CLAYTON) | (x$type == AC_CLAYTON)){
		.family = "C"
	}else
	if((x$type == HAC_FRANK) | (x$type == AC_FRANK)){
		.family = "F"
	}else
	if((x$type == HAC_JOE) | (x$type == AC_JOE)){
		.family = "J"
	}else
	if((x$type == HAC_AMH) | (x$type == AC_AMH)){
		.family = "A"
	}
	onacopulaL(.family, .tree2nacList(x$tree, .names, 1:length(.names)))
}

#------------------------------------------------------------------------------------------------------------------------

.tree2nacList = function(tree, names, numbers){	
	 if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
	 select = which(names %in% tree[-n])
	 
     if(any(s)){
         if(any(!s)){	
           res = list(tree[[n]], numbers[select], lapply(X = tree[which(!s)], FUN = .tree2nacList, names = names, numbers = numbers))
         }else{
           res = list(tree[[n]], numbers[select])
         }}else{
           res = list(tree[[n]], NULL, lapply(tree[which(!s)], FUN = .tree2nacList, names = names, numbers = numbers))
         }
     res
}