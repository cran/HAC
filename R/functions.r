# functions.r ############################################################################################################
# FUNCTION:         DESCRIPTION:
#  definitions      Global paramters.
#  TAU				Input argument of estimate.copula.
#  ML				Input argument of estimate.copula.
#  FML              Possible input argument of estimate.copula.
#  RML              Possible input argument of estimate.copula.
#  tau2theta       	Transforms Kendall's rank correlation coefficient to the dependence parameter of an Archimedean copula.
#  theta2tau        Transfrorms the dependence parameter of an Archimedean copula to Kendall's rank correlation coefficient.
#  phi        		The generator function of Archimedean copula.
#  phi.inv			The inverse of the generator function.
#  copMult        	Computes the value of d-dimensional AC for a given sample with values in [0,1]^d.
#  par.pairs        Returns the pairwise arranged parameter in a matrix, so that the parameters correspond to the lowest hierarchical level at which the variables are joined. 
#  .pair.matr       Supplementary function of par.pairs. Arranges the variables pairwise and returns the corresponding value. (Internal function)
##########################################################################################################################

 HAC_GUMBEL = 1
  AC_GUMBEL = 2
HAC_CLAYTON = 3
 AC_CLAYTON = 4
HAC_FRANK   = 5
 AC_FRANK   = 6
HAC_JOE     = 7
 AC_JOE     = 8
HAC_AMH     = 9
 AC_AMH     = 10

 ML = 1
FML = 2
RML = 3

#-------------------------------------------------------------------------------------------------------------------------------

tau2theta = function(tau, type = HAC_GUMBEL){
	n = length(tau)
	for(i in 1 : n){if((tau[i] < 0) | (tau[i] > 1)){return(warning(paste("tau[", i,"] should be in [0, 1).")))}}
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL))
        copGumbel@iTau(tau)
    else if((type == AC_CLAYTON) || (type == HAC_CLAYTON))
        copClayton@iTau(tau)
	else if((type == AC_FRANK) || (type == HAC_FRANK))
        copFrank@iTau(tau)
    else if((type == AC_JOE) || (type == HAC_JOE))
        copJoe@iTau(tau)
    else if((type == AC_AMH) || (type == HAC_AMH))
        copAMH@iTau(tau)
}

#-------------------------------------------------------------------------------------------------------------------------------

theta2tau = function(theta, type = HAC_GUMBEL){
	n = length(theta)
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		for(i in 1 : n){if(theta[i] < 1){return(warning(paste("theta[", i,"] >= 1 is required.")))}}
        copGumbel@tau(theta)
    }else if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
	    for(i in 1 : n){if(theta[i] <= 0){return(warning(paste("theta[", i,"] > 0 is required.")))}}
        copClayton@tau(theta)
    }else if((type == AC_FRANK) || (type == HAC_FRANK)){
	    for(i in 1 : n){if(theta[i] <= 0){return(warning(paste("theta[", i,"] > 0 is required.")))}}
        copFrank@tau(theta)
    }else if((type == AC_JOE) || (type == HAC_JOE)){
	    for(i in 1 : n){if(theta[i] < 1){return(warning(paste("theta[", i,"] >= 1 is required.")))}}
        copJoe@tau(theta)
    }else if((type == AC_AMH) || (type == HAC_AMH)){
	    for(i in 1 : n){if((theta[i] >= 1) || (theta[i] < 0)){return(warning(paste("theta[", i,"] >= 0 and < 1 is required.")))}}
        copAMH@tau(theta)
   }
}

#-------------------------------------------------------------------------------------------------------------------------------

phi = function(x, theta, type = HAC_GUMBEL){
	n = length(x)
	for(i in 1:n){if(x[i] < 0){return(warning(paste("x[", i,"] >= 0 is required.")))}}
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		if(theta >= 1){
			exp(-x^(1 / theta))
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		if(theta > 0){
			(x + 1)^(-1 / theta)
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == AC_FRANK) || (type == HAC_FRANK)){
		if(theta > 0){
			-log(1 - (1 - exp(-theta)) * exp(-x)) / theta
		}else{
			return(warning(paste("theta > 0 is required.")))
		}			
	}else
	if((type == AC_JOE) || (type == HAC_JOE)){
		if(theta >= 1){
			1 - (1 - exp(-x))^(1 / theta)
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}			
	}else
	if((type == AC_AMH) || (type == HAC_AMH)){
		if((theta >= 0) && (theta < 1)){
			(1 - theta) / (exp(x) - theta)
		}else{
			return(warning(paste("theta >= 0 and < 1 is required.")))
		}			
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

phi.inv = function(x, theta, type = HAC_GUMBEL){
	n = length(x)
	for(i in 1 : n){if((x[i] < 0) | (x[i] > 1)){return(warning(paste("x[", i,"] >= 0 and =< 1 is required.")))}}
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL)){
		if(theta >= 1){
			(-log(x))^theta
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		if(theta > 0){
			(x^(-theta) - 1)
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == AC_FRANK) || (type == HAC_FRANK)){
		if(theta > 0){
			-log((1 - exp(-x * theta))/(1 - exp(-theta)))
		}else{
			return(warning(paste("theta > 0 is required.")))
		}
	}else
	if((type == AC_JOE) || (type == HAC_JOE)){
		if(theta >= 1){
			-log(1 - (1 - x)^theta)
		}else{
			return(warning(paste("theta >= 1 is required.")))
		}
	}else
	if((type == AC_AMH) || (type == HAC_AMH)){
		if((theta >= 0) && (theta < 1)){
			log((1 - theta)/x + theta)
		}else{
			return(warning(paste("theta >= 0 and < 1 is required.")))
		}			
	}
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
           if(FUN == "TAU")
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
     
     if(any(s)){
         if(any(!s)){
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

