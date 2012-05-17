# copula_functions.r #####################################################################################################
# FUNCTION:         DESCRIPTION:
#  definitions      Global paramters.
#  theta.eps		Is addressed during the estimation procedure.
#  TAU				Possible input argument of estimate.copula.
#  ML				Possible input argument of estimate.copula.
#  FML              Possible input argument of estimate.copula.
#  tau2theta       	Convertes Kendall's rank correlation coefficient into dependence parameter.
#  theta2tau        Convertes dependence parameter into Kendall's rank correlation coefficient.
#  phi        		The generator function.
#  phi.inv			The inverse of the generator function.
#  copMult        	Computes the value of d-dimensional AC.
#  par.pairs        Returns the pairwise arranged parameter in a matrix.
#  .pair.matr       Supplementary function of par.pairs. (Internal function)       
##########################################################################################################################

HAC_GUMBEL = 0
AC_GUMBEL = 1
HAC_CLAYTON = 3
AC_CLAYTON = 4
GAUSS = 5

TAU = 0
ML  = 1
FML = 2

#-------------------------------------------------------------------------------------------------------------------------------

tau2theta = function(tau, type = HAC_GUMBEL){
	n = length(tau)
	if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		for(i in 1 : n){if((tau[i] < 0) | (tau[i] > 1)){return(warning(paste("Element[", i,"] should be in [0, 1).")))}}
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL))
        1 / (1-tau)
    else if((type == AC_CLAYTON) || (type == HAC_CLAYTON))
        2 * tau/ (1-tau)}
    else if(type == GAUSS){
    	for(i in 1 : n){if(abs(tau[i]) > 1){return(warning(paste("Element[", i,"] should be in [-1, 1].")))}}
        sin(tau * pi / 2)}
}

#-------------------------------------------------------------------------------------------------------------------------------

theta2tau = function(theta, type = HAC_GUMBEL){
	n = length(theta)
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		for(i in 1 : n){if(theta[i] < 1){return(warning(paste("Element[", i,"] >= 1 is required.")))}}
        1 - 1/theta}
    else if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
	    for(i in 1 : n){if(theta[i] <= 0){return(warning(paste("Element[", i,"] > 0 is required.")))}}
        theta / (2 + theta)}
    else if(type == GAUSS){
	    for(i in 1 : n){if(abs(theta[i]) > 1){return(warning(paste("Element[", i,"] < |1| is required.")))}}
        2/pi * asin(theta)}
}

#-------------------------------------------------------------------------------------------------------------------------------

phi = function(x, theta, type = HAC_GUMBEL){
	n = length(x)
	for(i in 1:n){if(x[i] < 0){return(warning(paste("Element[", i,"] >= 0 is required.")))}}
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		if(theta >= 1){
			exp(-x^(1 / theta))}
		else
			return(warning(paste("theta >=", 1," is required.")))}
	else{
	if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		if(theta > 0){
			(x + 1)^(-1 / theta)}
		else
		return(warning(paste("theta >=", 0," is required.")))}}
}

#-------------------------------------------------------------------------------------------------------------------------------

phi.inv = function(x, theta, type = HAC_GUMBEL){
	n = length(x)
	for(i in 1 : n){if((x[i] < 0) | (x[i] > 1)){return(warning(paste("Element[", i,"] >= 0 and =< 1 is required.")))}}
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		if(theta >= 1){
			(-log(x))^theta}
		else
			return(warning(paste("theta >=", 1," is required.")))}
	else{
	if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		if(theta > 0){
			(x^(-theta) - 1)}
		else
			return(warning(paste("theta >", 0," is required.")))}}
}

#-------------------------------------------------------------------------------------------------------------------------------

copMult = function(X, theta = 1, type = HAC_GUMBEL){	
	phi(rowSums(phi.inv(X, theta = theta, type)), theta = theta, type)
}

#------------------------------------------------------------------------------------------------------------------------

par.pairs = function(hac, FUN = NULL, ...){
    tree = hac$tree
    vars = .get.leaves(tree)
    d = length(vars)
    matr = matrix(NA,nrow=d,ncol=d); colnames(matr)=rownames(matr)=vars
    matr = .pairs.matr(tree, matr)
    diag(matr) = NA
    
    if(class(FUN)=="function"){
        matr[lower.tri(matr)] = FUN(matr[lower.tri(matr)], ...)
        matr[upper.tri(matr)] = FUN(matr[upper.tri(matr)], ...)
    }else{
        if(!is.null(FUN)){
           if(FUN==TAU)
            matr[lower.tri(matr)] = theta2tau(matr[lower.tri(matr)], type=hac$type)
            matr[upper.tri(matr)] = theta2tau(matr[upper.tri(matr)], type=hac$type)
        }}
    matr
}

#------------------------------------------------------------------------------------------------------------------------

.pairs.matr = function(tree, matr){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     
     if(any(s==TRUE)){
         if(any(s==FALSE)){
            l = sapply(tree[which(!s)], .get.leaves)
            for(i in 1:(length(l)-1))for(j in (i+1):length(l))matr[unlist(l[[i]]),unlist(l[[j]])]=matr[unlist(l[[j]]),unlist(l[[i]])]=tree[[n]]
            for(i in 1:length(l))matr[unlist(l[[i]]),unlist(tree[which(s)])]=matr[unlist(tree[which(s)]),unlist(l[[i]])]=tree[[n]]
            matr[unlist(tree[which(s)]), unlist(tree[which(s)])]=tree[[n]]
            for(i in which(!s)){
                matr = .pairs.matr(tree[i], matr)         
            }
         }else{
            matr[unlist(tree[-n]), unlist(tree[-n])]=tree[[n]]
         }}else{
            l = sapply(tree[-n], .get.leaves)
            for(i in 1:(length(l)-1))for(j in (i+1):length(l))matr[unlist(l[[i]]),unlist(l[[j]])]=matr[unlist(l[[j]]),unlist(l[[i]])]=tree[[n]]
            for(i in 1:(n-1)){
                 matr = .pairs.matr(tree[i], matr)            
            }
        }
    return(matr)
}
