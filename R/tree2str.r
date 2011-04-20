# tree2str.r #############################################################################################################
# FUNCTION:               	DESCRIPTION:
#  tree2str					Prints the structure of HAC as string.
#  .tree2str				Returns the structure of nested Archimedean copulae. (Internal function)
#  get.params				Prints the dependency parameters.
#  .get.params				Returns the dependency parameter of nested Archimedean copulae. (Internal function)
##########################################################################################################################

tree2str = function(hac, theta = TRUE, digits = 3){
	
    if((hac$type == HAC_GUMBEL) || (hac$type == HAC_ROTATED_GUMBEL) || (hac$type == HAC_CLAYTON)){
        res = .tree2str(hac$model, theta, digits)
    }else if((hac$type == AC_GUMBEL) || (hac$type == AC_CLAYTON)){
        res = paste("(", paste(paste("V", 1:hac$model$dim, sep = ""), collapse = "."), ")", 
        if(theta){paste("_", round(hac$model$theta, digits = digits), sep = "")}else{""}, sep = "")
    }else if(hac$type == GAUSS){
        res = "No tree for GAUSS models"
    }
    res
}

#-------------------------------------------------------------------------------------------------------------------------------

.tree2str = function(tree, theta = TRUE, dec = 3){
	if((class(tree$V1) == "character") & (class(tree$V2) == "character")){
		if(theta){return(paste("(", paste(tree$V1, collapse = "."), ".", paste(tree$V2, collapse = "."), ")_{", round(tree$theta, digits = dec), "}", sep = ""))}
		else{return(paste("(", paste(tree$V1, collapse = "."), ".", paste(tree$V2, collapse = "."), ")", sep = ""))}
	}
	else if((class(tree$V1) != "character") & (class(tree$V2) == "character")){
		if(theta){return(paste("(", .tree2str(tree$V1, theta, dec), ".", paste(tree$V2, collapse = "."), ")_{", round(tree$theta, digits = dec), "}", sep = ""))}
		else{return(paste("(", .tree2str(tree$V1, theta, dec), ".", paste(tree$V2, collapse = "."), ")", sep = ""))}
	}
	else if((class(tree$V1) == "character") & (class(tree$V2) != "character")){
		if(theta){return(paste("(", paste(tree$V1, collapse = "."), ".", .tree2str(tree$V2, theta, dec), ")_{", round(tree$theta, digits = dec), "}", sep = ""))}
		else{return(paste("(", paste(tree$V1, collapse = "."), ".", .tree2str(tree$V2, theta, dec), ")", sep = ""))}
	}
	else if((class(tree$V1) != "character") & (class(tree$V2) != "character")){
		if(theta){return(paste("(", .tree2str(tree$V1, theta, dec), ".", .tree2str(tree$V2, theta, dec), ")_{", round(tree$theta, digits = dec), "}", sep = ""))}
		else{return(paste("(", .tree2str(tree$V1, theta, dec), ".", .tree2str(tree$V2, theta, dec), ")", sep = ""))}
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

get.params = function(hac){
	
    if((hac$type == HAC_GUMBEL) || (hac$type == HAC_ROTATED_GUMBEL) || (hac$type == HAC_CLAYTON)){
        res = .get.params(hac$model)
    }else if((hac$type == AC_GUMBEL) || (hac$type == AC_CLAYTON)){
        res = hac$model$theta
    }else if(hac$type == GAUSS){
        res = hac$model
    }
    res
}

#-------------------------------------------------------------------------------------------------------------------------------

.get.params = function(hac){
	if((class(hac$V1) == "character") & (class(hac$V2) == "character")){
		return(hac$theta)
	}
	else if((class(hac$V1) != "character") & (class(hac$V2) == "character")){
		return(c(.get.params(hac$V1), hac$theta))
	}
	else if((class(hac$V1) == "character") & (class(hac$V2) != "character")){
		return(c(hac$theta, .get.params(hac$V2)))
	}
	else if((class(hac$V1) != "character") & (class(hac$V2) != "character")){
		return(c(hac$theta, .get.params(hac$V1), .get.params(hac$V2)))
	}
}