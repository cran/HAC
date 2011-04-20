# copula_cdf.r ###########################################################################################################
# FUNCTION:               	DESCRIPTION:
#  pHAC						Computes the values of the cdf for a given sample and 'hac' object.
#  .cop.cdf           		A recursive function. Computes the values of the cdf for a given type of copula. (Internal function)
#  .get.names.in.tree       Returns the columns of the sample for which the cdf is to compute. (Internal function)           
##########################################################################################################################

pHAC = function(X, hac){
    if(hac$type == HAC_ROTATED_GUMBEL){
        return((1 - .cop.cdf(hac$model, X, HAC_GUMBEL)))
    }else if((hac$type == HAC_GUMBEL) || (hac$type == HAC_CLAYTON)){
        return(.cop.cdf(hac$model, X, hac$type))
    }else if(hac$type == GAUSS){
        return(pcopula(normalCopula(hac$model[lower.tri(hac$model)], dim = dim(X)[2], dispstr = "un"), X))
    }else if(hac$type == AC_GUMBEL){
        return(pcopula(gumbelCopula(hac$model$theta, dim = dim(X)[2]), X))
    }else if(hac$type == AC_CLAYTON){
        return(pcopula(claytonCopula(hac$model$theta, dim = dim(X)[2]), X))
    }
}

#------------------------------------------------------------------------------------------------------------------------

.cop.cdf = function(Ltree, sample, copula_type){
	if((class(Ltree$V1) == "character") & (class(Ltree$V2) == "character")){
		return(copMult(cbind(sample[,Ltree$V1], sample[,Ltree$V2]), Ltree$theta, copula_type))
	}else if((class(Ltree$V1) != "character") & (class(Ltree$V2) == "character")){
		return(copMult(cbind(.cop.cdf(Ltree$V1, sample[,.get.names.in.tree(Ltree$V1)], copula_type), sample[,Ltree$V2]), Ltree$theta, copula_type))
	}else if((class(Ltree$V1) == "character") & (class(Ltree$V2) != "character")){
		return(copMult(cbind(sample[,Ltree$V1], .cop.cdf(Ltree$V2, sample[,.get.names.in.tree(Ltree$V2)], copula_type)), Ltree$theta, copula_type))
	}else if((class(Ltree$V1) != "character") & (class(Ltree$V2) != "character")){
		return(copMult(cbind(.cop.cdf(Ltree$V1, sample[,.get.names.in.tree(Ltree$V1)], copula_type), .cop.cdf(Ltree$V2, sample[,.get.names.in.tree(Ltree$V2)], copula_type)), Ltree$theta, copula_type))
	}
}

#------------------------------------------------------------------------------------------------------------------------
	
.get.names.in.tree = function(Ltree){
	if((class(Ltree$V1) == "character") & (class(Ltree$V2) == "character")){
		return(c(Ltree$V1, Ltree$V2))
	}else if((class(Ltree$V1) != "character") & (class(Ltree$V2) == "character")){
		return(c(.get.names.in.tree(Ltree$V1), Ltree$V2))
	}else if((class(Ltree$V1) == "character") & (class(Ltree$V2) != "character")){
		return(c(Ltree$V1, .get.names.in.tree(Ltree$V2)))
	}else if((class(Ltree$V1) != "character") & (class(Ltree$V2) != "character")){
		return(c(.get.names.in.tree(Ltree$V1), .get.names.in.tree(Ltree$V2)))
	}
}
	