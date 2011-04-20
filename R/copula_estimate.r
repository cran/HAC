# copula_estimate.r ######################################################################################################
# FUNCTION:               	DESCRIPTION:
#  estimate.copula			Estimates the structure and the parameter of a HAC for a given sample X.
#  .ub          			Assures that the dependency parameter of the initial node is smaller than parameter of consecutive nodes. (Internal function)       
##########################################################################################################################

estimate.copula = function(X, type = HAC_GUMBEL, method = TAU, epsilon = 0){
	
	if(is.null(colnames(X)) == TRUE){colnames(X) = paste("X", 1 : dim(X)[2], sep = "")}
	
    if((method == TAU) && (dim(X)[1] == 1))stop("Cant estimate copula based on the tau method with just one observation")
    if((type == HAC_GUMBEL) || (type == HAC_CLAYTON)){
        main.dim = dim(X)[2]
        tree = as.list(colnames(X))
    	for(main.i in 1:(main.dim-2)){
    
            if(method == TAU){
            X.help = cor(X, method = "kendall")
            X.help[which(X.help < 0)] = 0
            matr = tau2theta(X.help, type)
            }else{
                matr = matrix(0, (main.dim-main.i+1), (main.dim-main.i+1))
                for(i in 1:((main.dim-main.i+1)-1))for(j in (i+1):(main.dim-main.i+1))
                    matr[i, j] = matr[j, i] = tau2theta(optimise(f = function(y, i, j){sum(log(dAC(X[,i], X[,j], tau2theta(y, type), type)))}, i = i, j = j, interval = c(0 + theta.eps, 1 - theta.eps), maximum = TRUE)$maximum, type)
            }

            cur.dim = dim(matr)[1]
            max.m = -10
            max.i = max.ii = 0
            for(i in 1:cur.dim)for(ii in i:cur.dim)if(i != ii){if(matr[i,ii] > max.m){max.m = matr[i,ii];max.i = i;max.ii = ii}} # what should be joined

            sub.min = 1000.
            if(class(tree[[max.i]]) != "character") sub.min = c(sub.min, tree[[max.i]]$theta)
            if(class(tree[[max.ii]]) != "character") sub.min = c(sub.min, tree[[max.ii]]$theta)
            if(min(min(sub.min), matr[max.i,max.ii]) == matr[max.i,max.ii]){sub.min = matr[max.i,max.ii]}else{sub.min = min(sub.min) - theta.eps}

            co = copMult(cbind(X[,max.i], X[,max.ii]), max(sub.min, tau2theta(theta.eps)), type)
            X = matrix(X[,-max(max.i, max.ii)], ncol = (main.dim-main.i))
            X[,min(max.i, max.ii)] = co

            tree[[max.i]] = list(V1 = tree[[max.i]], V2 = tree[[max.ii]], theta = max(sub.min, tau2theta(theta.eps)))
            tree = tree[-max.ii]
        }
        if(method == TAU){
            res = list(V1 = tree[1][[1]], V2 = tree[2][[1]], theta = tau2theta(max(cor(X[,1], X[,2], method = "kendall"), theta.eps), type))
            res = .union(res, epsilon)
        }else{  			
            res = list(V1 = tree[1][[1]], V2 = tree[2][[1]], theta = tau2theta(optimise(f = function(y, i, j){sum(log(dAC(X[,1], X[,2], tau2theta(y, type), type)))}, i = 1, j = 2, interval = c(0 + theta.eps, .ub(tree[1][[1]], tree[2][[1]], type)), maximum = TRUE)$maximum, type))
            res = .union(res, epsilon)
        }
    }else if(type == GAUSS){
        res = cov(qnorm(X))
        return(res)
    }else if(type == AC_GUMBEL){
        if(method == TAU){
            stop("No estimation based on taus for simple Archimedean copulas (with d>2, (d=2 not implemented yet))")
        }else{
            res = list(dim = dim(X)[2], theta = fitCopula(gumbelCopula(1.5, dim = dim(X)[2]), X, method = "ml")@estimate)
        }
    }else if(type == AC_CLAYTON){
        if(method == TAU){
            stop("No estimation based on taus for simple Archimedean copulas (with d>2, (d=2 not implemented yet))")
        }else{
            res = list(dim = dim(X)[2], theta = fitCopula(claytonCopula(1.5, dim = dim(X)[2]), X, method = "ml")@estimate)
        }
    }else if(type == HAC_ROTATED_GUMBEL){
        stop("no yet implemented")
    }
    result = list(type = type, model = res)
	class(result) = "hac"
	return(result)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.ub = function(tree.1, tree.2, type){
  	if((class(tree.1) != "character") & (class(tree.2) == "character"))
  		theta2tau(tree.1$theta, type)
  	else 
  	if((class(tree.1) == "character") & (class(tree.2) != "character"))
  		theta2tau(tree.2$theta, type)
	else 
	if((class(tree.1) != "character") & (class(tree.2) != "character"))
  		theta2tau(min(tree.1$theta, tree.2$theta), type)
}