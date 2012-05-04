# copula_estimate.r ######################################################################################################
# FUNCTION:               	DESCRIPTION:
#  estimate.copula			Estimates the structure and the parameter of a HAC for a given sample.
#  .ub         			 	Assures the dependency parameter of the initial node being smaller than parameter of consecutive nodes. (Internal function)
#  . margins				Estimates the d marginal distributions for a d dimensional sample. (Internal function)   
#  . one.mar				Estimates one marginal distributions for a given sample. (Internal function)   
#  .max.min					0's contained in the data matrix are set to 0.000001 and 1's to 1-000001. (Internal function)
#  .constraints.ui          Returns the matrix for the constraints, i.e. ui of constrOptim.
#  .constraints.ui          Returns the vector of lower bounds, i.e. ci of constrOptim.
#  .rebuild                 Creates the tree of an hac-object for an ordered parameter vector.
##########################################################################################################################

estimate.copula = function(X, type = HAC_GUMBEL, method = ML, hac = NULL, epsilon = 0, agg.method = "mean", margins = NULL, theta.eps = 0.001, na.rm = FALSE, max.min = TRUE, ...){
	
	if(is.null(colnames(X))){g.names = names = paste("X", 1 : NCOL(X), sep = "")}else g.names = names = colnames(X)
	
	X = .margins(X, margins)
	colnames(X) = names
	
	if(na.rm == TRUE){
		X =  na.omit(X, ...)}

	if(max.min == TRUE){
		X = .max.min(X)}
		
    if((method == TAU) && (dim(X)[1] == 1))stop("Cant estimate copula based on the tau method with just one observation")
    if(((type == HAC_GUMBEL) | (type == HAC_CLAYTON)) & (NCOL(X)>2)){
    
        if(method != FML){
           main.dim = NCOL(X); tree = as.list(names)
    	
        for(main.i in 1:(main.dim-2)){
    
            if(method == TAU){
            X.help = cor(X, method = "kendall")
            X.help[which(X.help < 0)] = theta.eps
            matr = tau2theta(X.help, type)
            }else{
                matr = matrix(0, (main.dim-main.i+1), (main.dim-main.i+1))
                for(i in 1:((main.dim-main.i+1)-1))for(j in (i+1):(main.dim-main.i+1))
                    matr[i, j] = matr[j, i] = tau2theta(optimise(f = function(y, i, j){sum(log(.dAC(X[,i], X[,j], tau2theta(y, type), type)))}, i = i, j = j, interval = c(0 + theta.eps, 1 - theta.eps), maximum = TRUE)$maximum, type)
            }

            cur.dim = NROW(matr); max.m = -10;max.i = max.ii = 0; sub.min = 1000
            for(i in 1:cur.dim)for(ii in i:cur.dim)if(i != ii){if(matr[i,ii] > max.m){max.m = matr[i,ii];max.i = i;max.ii = ii}}

            if(class(tree[[max.i]]) != "character") sub.min = c(sub.min, tree[[max.i]][[length(tree[[max.i]])]])
            if(class(tree[[max.ii]]) != "character") sub.min = c(sub.min, tree[[max.ii]][[length(tree[[max.ii]])]])
            if(min(min(sub.min), matr[max.i,max.ii]) == matr[max.i,max.ii]){sub.min = matr[max.i,max.ii]}else{sub.min = min(sub.min) - theta.eps}

            co = copMult(cbind(X[,max.i], X[,max.ii]), max(sub.min, tau2theta(theta.eps, type)), type)
            X = matrix(X[,-max(max.i, max.ii)], ncol = (main.dim-main.i))
            X[,min(max.i, max.ii)] = co

            tree[[max.i]] = list(tree[[max.i]], tree[[max.ii]], max(sub.min, tau2theta(theta.eps, type)))
            tree = tree[-max.ii]; main.i = main.i+1
            }
        
        if(method == TAU){
            res = c(list(tree[[1]]), list(tree[[2]]), tau2theta(max(cor(X[,1], X[,2], method = "kendall"), theta.eps), type))
            res = .union(res, epsilon, method = agg.method, ...)
        }else{  			
            res = c(list(tree[[1]]), list(tree[[2]]), tau2theta(optimise(f = function(y, i, j){sum(log(.dAC(X[,1], X[,2], tau2theta(y, type), type)))}, i = 1, j = 2, interval = c(0 + theta.eps, .ub(tree[[1]][[length(tree[[1]])]], tree[[2]][[length(tree[[2]])]], type)), maximum = TRUE)$maximum, type))
            res = .union(res, epsilon, method = agg.method, ...)
    }}else{
       if(is.null(hac)){
          stop("A hac object is required.")
       }else{
          values = get.params(hac, sort.v = TRUE, decreasing=FALSE)
          tree.full = hac$tree
          n.params = length(values)
          initial = numeric(n.params); initial[1]=if(type == HAC_GUMBEL){1}else{0}
          constr = .constraints.ui(tree.full)
          m = diag(1, nrow = constr+1, ncol = constr+1)
          LL = to.logLik(X, hac)
          optim = constrOptim(values, f=LL, grad=NULL, ui=m, ci=as.vector(.constraints.ci(tree.full, values, init = initial)), control=list(fnscale=-1), hessian=FALSE)
          theta = optim$par
          res = .rebuild(tree.full, values, theta)
    }}}else{ 
    
    if((type == AC_GUMBEL) | (type == HAC_GUMBEL)){
        if(method == TAU){
            stop("No estimation based on taus for simple Archimedean copulas (with d>2, (d=2 not implemented yet))")
        }else{
            res = c(as.list(names), fitCopula(gumbelCopula(1.5, dim = dim(X)[2]), X, method = "ml")@estimate)
        }
    }else if((type == AC_CLAYTON) | (type == HAC_CLAYTON)){
        if(method == TAU){
            stop("No estimation based on taus for simple Archimedean copulas (with d>2, (d=2 not implemented yet))")
        }else{
            res = c(as.list(names), fitCopula(claytonCopula(1.5, dim = dim(X)[2]), X, method = "ml")@estimate)
        }
    }else if(type == HAC_ROTATED_GUMBEL){
        stop("no yet implemented")
    }}
	hac(type = type, tree = res)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.ub = function(tree.1, tree.2, type){
  	if((class(tree.1) == "numeric") & (class(tree.2) == "character"))
  		theta2tau(tree.1, type)
  	else 
  	if((class(tree.1) == "character") & (class(tree.2) == "numeric"))
  		theta2tau(tree.2, type)
	else 
	if((class(tree.1) == "numeric") & (class(tree.2) == "numeric"))
  		theta2tau(min(tree.1, tree.2), type)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.margins = function(X, margins, ...){
	if(is.null(margins)){
		X
	}else{
		if(length(margins)==1){
		X = apply(X, 2, .one.mar, spec = margins, ...)
	}else{
		for(i in 1:NCOL(X)){X[,i] = .one.mar(X[,i], margins[i],...)}}
		X
	}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.one.mar = function(X, spec, ...){
	n = NROW(X)
		if(spec == "edf"){
			f = ecdf(X, ...)
			n/(n+1)*f(X)
		}
		#else{
		#if(spec == "np"){
		#	fitted(np::npudist(~X, ...))
		#}
	else{
		.opt.margin(data = X, spec = spec)}
	#}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.opt.margin = function(data, spec){
	boundary = 10000
	if((spec == "beta") | (spec == "cauchy") | (spec == "chisq") | (spec == "f") | (spec == "gamma") | (spec == "lnorm") | (spec == "norm") | (spec == "t") | (spec == "weibull")){
		loglik = function(par, data){sum(log(eval(do.call(paste("d", spec, sep = ""), args = list(x = data, par[1], par[2])))))}
		op =constrOptim(theta = c(1, 1), f = loglik, grad = NULL, ui = matrix(c(1, 0, -1, 0, 0, 1), nrow = 3, byrow = TRUE), ci = c(-rep(boundary, 2), 0), data = data, control = list(fnscale = -1), hessian = FALSE)
		eval(do.call(paste("p", spec, sep = ""), args = list(q = data, op$par[1], op$par[2])))					
	}else{
	if((spec == "exp")){
		op =optimise(f = function(par, data){sum(log(eval(do.call(paste("d", spec, sep = ""), args = list(x = data, par)))))}, data = data
, lower = 0.0001, upper = 100, maximum = TRUE)$maximum
		eval(do.call(paste("p", spec, sep = ""), args = list(q = data, op)))	
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.max.min = function(X){
		if((any(X<=0)) |  (any(X>=1))){
			X[which(X >= 1)] = 0.999999
			X[which(X <= 0)] = 0.000001
			X}
		else{
			X}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.constraints.ui = function(tree){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
 
     if(any(s==TRUE)){
         if(any(s==FALSE)){
            n.constr = length(which(!s))
            n.constr = sum(n.constr, sapply(tree[which(!s)], .constraints.ui))
         }else{
            n.constr = 0
         }}else{
            n.constr = sum(length(which(!s)), sapply(tree[-n], .constraints.ui))
        }
    return(n.constr)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.constraints.ci = function(tree, values, init){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
 
     if(any(s==TRUE)){
         if(any(s==FALSE)){
            n.constr = length(which(!s))
            init[which(values %in% sort(sapply(tree[which(!s)], FUN=function(r){r[[length(r)]]})))] = rep(tree[[n]], n.constr)
            
            for(i in which(!s)){
                init = .constraints.ci(tree[i], values, init)            
            }
         }}else{
            n.constr = length(s)
            init[which(values %in% sort(sapply(tree[-n], FUN=function(r){r[[length(r)]]})))] = rep(tree[[n]], n.constr)
            
            for(i in 1:(n-1)){
                 init = .constraints.ci(tree[i], values, init)            
            }
        }
    return(init)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.rebuild = function(tree, values, theta){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     tree[[n]]=theta[which(values==tree[[n]])]
                 
     if(any(s==TRUE)){
         if(any(s==FALSE)){
            tree=c(tree[which(s)], lapply(tree[which(!s)], .rebuild, values, theta), tree[[n]])           
        }else{
            tree=tree
        }}else{
            tree = c(lapply(tree[-n], .rebuild, values, theta), tree[[n]])
        }
    return(tree)
}
