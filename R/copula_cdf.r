# copula_cdf.r ###########################################################################################################
# FUNCTION: 			DESCRIPTION:
#  pHAC					Computes the values of the cdf for a given sample and 'hac' object.
#  .cop.cdf           	Supplementary function for pHAC. (Internal function)
##########################################################################################################################

pHAC = function(X, hac, margins = NULL, na.rm = FALSE, ...){
	
	X = .margins(X, margins)
	
	if(na.rm == TRUE){
			X = na.omit(X, ...)}
			
    if(hac$type == HAC_ROTATED_GUMBEL){
        return((1 - .cop.cdf(X, hac$tree, HAC_GUMBEL)))
    }else if((hac$type == HAC_GUMBEL) || (hac$type == HAC_CLAYTON)){
        return(.cop.cdf(X, hac$tree, hac$type))
    }else if(hac$type == GAUSS){
        return(pcopula(normalCopula(hac$tree[lower.tri(hac$tree)], dim = NCOL(X), dispstr = "un"), X))
    }else if(hac$type == AC_GUMBEL){
    	n = length(unlist(hac$tree))
        return(pcopula(gumbelCopula(hac$tree[[n]], dim = (n-1)), X))
    }else if(hac$type == AC_CLAYTON){
    	n = length(unlist(hac$tree))
        return(pcopula(claytonCopula(hac$tree[[n]], dim = (n-1)), X))
    }
}

#------------------------------------------------------------------------------------------------------------------------

.cop.cdf = function(sample, tree, type){
	if(length(tree)==1){tree = tree[[1]]}
	n = length(tree); m = matrix(0, nrow = NROW(sample)); names = colnames(sample)
	s = sapply(tree, is.character)

	if(any(s[-n]==TRUE)){
		if(any(s[-n]==FALSE)){
			select = unlist(tree[s])
				for(i in 1:length(select)){select[i]=(which(names==select[i]))}; select = as.numeric(select)
				exclude = c(1:(n-1))[which(s[-n]==FALSE)]
			m = copMult(cbind(sample[, select], sapply(tree[exclude], .cop.cdf, sample = sample[, -select], type = type)), theta = tree[[n]], type = type)
		}else{
			m = copMult(cbind(sample[, unlist(tree[s])]), theta = tree[[n]], type = type)
	}}else{
		m = copMult(sapply(tree[-n], .cop.cdf, sample = sample, type = type), theta = tree[[n]], type = type)
	}
	return(m)
}
