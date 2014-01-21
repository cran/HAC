# estimate.r ############################################################################################################# 
# FUNCTION:               	DESCRIPTION: 
#  estimate.copula			Estimates the structure and the parameter of a HAC for a given sample. 
#  .ML.TAU                  Estimation procedures based on binary trees, i.e., for method = ML and method = TAU. (Internal function)
#  .FML                     Full Maximum Likelihood (FML) estimation procedure. It needs an 'hac' object as argument to construct the log-likelihood which depends on the structure of the HAC. (Internal function)
#  .RML                     Recursive Maximum Likelihood (RML) estimation procedure as discussed in Okhrin, Okhrin and Schmid (2013). (Internal function)
#  .ub         			 	Enures the dependency parameter of the initial node being smaller than parameter of consecutive nodes. (Internal function) 
#  .margins				    Estimates the marginal distributions and returns the fitted values for a d-dimensional sample. (Internal function)   
#  .one.mar				    Estimates one marginal distributions for a given univariate sample. (Internal function)   
#  .max.min					0's contained in the data matrix are set to 0.000001 and 1's to 1-0.00001. (Internal function) 
#  .constraints.ui          Returns a matrix of constraints according to the matrix ui of constrOptim. This matrix ensures the parameters being increasing from the highest to the lowest hierarchical level for the full ML approach. (Internal function)
#  .rebuild                 Matches the tree of a 'hac' object according to an ordered parameter vector. (Internal function) 
########################################################################################################################## 

estimate.copula = function(X, type = HAC_GUMBEL, method = RML, hac = NULL, epsilon = 0, agg.method = "mean", margins = NULL, na.rm = FALSE, max.min = TRUE, ...){
	
	if(is.null(colnames(X))){g.names = names = paste("X", 1 : NCOL(X), sep = "")}else names = colnames(X)
	
	X = .margins(X, margins)
	colnames(X) = names
	
	if(na.rm){X = na.omit(X, ...)}

	if(max.min){X = .max.min(X)}
	
	d = NCOL(X)	
    if(((type == HAC_GUMBEL) | (type == HAC_CLAYTON) | (type == HAC_FRANK) | (type == HAC_JOE) | (type == HAC_AMH)) & (d > 2)){
    
        if(method == ML){
            res = .ML(X = X, type = type, epsilon = epsilon, agg.method = agg.method, names = names, ...)    
        }else{
        if(method == FML){
            if(is.null(hac)){
                stop("A hac object is required.")
            }else{
                res = .FML(X = X, type = type, hac = hac)
        	}
        }else{
            res = .RML(X = X, type = type, epsilon = epsilon, agg.method = agg.method, names = names, ...)
        }}
    }else{
    	if((type == AC_GUMBEL) | (type == HAC_GUMBEL)){
            res = c(as.list(names), fitCopula(gumbelCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
    	}else 
    	if((type == AC_CLAYTON) | (type == HAC_CLAYTON)){
            res = c(as.list(names), fitCopula(claytonCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == AC_FRANK) | (type == HAC_FRANK)){
            res = c(as.list(names), fitCopula(frankCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == AC_JOE) | (type == HAC_JOE)){
            res = c(as.list(names), fitCopula(joeCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == AC_AMH) | (type == HAC_AMH)){
            res = c(as.list(names), fitCopula(amhCopula(param = 0.25, dim = d), X, method = "ml")@estimate)
        }
    }
	hac(type = type, tree = res)
}     

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.ML = function(X, type, epsilon = 0, agg.method = "mean", names, ...){
        main.dim = NCOL(X); tree = as.list(names); upper.b = if((type == AC_AMH) | (type == HAC_AMH)){1/3-1e-8}else{1-1e-8}
        for(main.i in 1:(main.dim-2)){
    
           matr = matrix(0, (main.dim-main.i+1), (main.dim-main.i+1))
                for(i in 1:((main.dim-main.i+1)-1))for(j in (i+1):(main.dim-main.i+1))
                    matr[i, j] = matr[j, i] = tau2theta(optimise(f = function(y, i, j){sum(log(.dAC(X[,i], X[,j], tau2theta(y, type), type)))}, i = i, j = j, interval = c(1e-8, upper.b), maximum = TRUE)$maximum, type)
            

            cur.dim = NROW(matr); max.m = -10; max.i = max.ii = 0; sub.min = 1000
            for(i in 1:cur.dim)for(ii in i:cur.dim)if(i != ii){if(matr[i,ii] > max.m){max.m = matr[i,ii];max.i = i;max.ii = ii}}

            if(class(tree[[max.i]]) != "character") sub.min = c(sub.min, tree[[max.i]][[length(tree[[max.i]])]])
            if(class(tree[[max.ii]]) != "character") sub.min = c(sub.min, tree[[max.ii]][[length(tree[[max.ii]])]])
            if(min(min(sub.min), matr[max.i,max.ii]) == matr[max.i,max.ii]){sub.min = matr[max.i,max.ii]}else{sub.min = min(sub.min) - 1e-8}

            co = copMult(cbind(X[,max.i], X[,max.ii]), max(sub.min, tau2theta(1e-8, type)), type)

            X = matrix(X[,-max(max.i, max.ii)], ncol = (main.dim-main.i))
            X[,min(max.i, max.ii)] = co

            tree[[max.i]] = list(tree[[max.i]], tree[[max.ii]], max(sub.min, tau2theta(1e-8, type)))
            tree = tree[-max.ii]; main.i = main.i+1
            }
         			
            res = c(list(tree[[1]]), list(tree[[2]]), tau2theta(optimise(f = function(y){sum(log(.dAC(X[,1], X[,2], tau2theta(y, type), type)))}, interval = c(1e-8, .ub(tree[[1]][[length(tree[[1]])]], tree[[2]][[length(tree[[2]])]], type)), maximum = TRUE)$maximum, type))
    
    .union(res, epsilon = epsilon, method = agg.method, ...)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.FML = function(X, type, hac){
    values = get.params(hac, sort.v = TRUE, decreasing=FALSE)
    tree.full = hac$tree
	initial=if((type == HAC_GUMBEL) | (type == HAC_JOE)){1+1e-8}else{1e-8}
    ui = .constraints.ui(tree.full, m = matrix(c(1, rep(0, length(values)-1)), nrow=1), values = values)
    LL = to.logLik(X, hac)
    optim = constrOptim(values, f=LL, grad=NULL, ui=as.matrix(ui), ci=as.vector(c(initial, rep(1e-8, NROW(ui)-1))), control=list(fnscale=-1), hessian=FALSE)
    .rebuild(tree.full, values, optim$par)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.RML = function(X, type, epsilon, agg.method, names, ...){
    main.dim = NCOL(X); tree = as.list(names)
    select = if((type == AC_AMH) | (type == HAC_AMH)){c(0, 0, 1-1e-8)}else{c(0, 0, 100)}
        
        while(main.dim > 1){       
           matr.p = matrix(0, main.dim, main.dim)
           ff.done = NULL
                for(i in 1:(main.dim-1)){
                   for(j in (i+1):main.dim){
                        if((names[i] != "tree") & (names[j] != "tree")){
                            matr.p[i, j] = matr.p[j, i] = optimise(f = function(y, i, j){sum(log(.dAC(X[, names[i]], X[, names[j]], tau2theta(y, type), type)))}, i = i, j = j, interval = c(1e-8, theta2tau(select[3], type)), maximum = TRUE)$maximum
                        }else{
                        if(((names[i] == "tree") & (names[j] != "tree")) | ((names[i] != "tree") & (names[j] == "tree"))){
                            for(l in 1:length(tree[-without])){
                                if(is.null(ff.done)){
                                    tree.nn = list("leaf.new", .tree.without.params(tree[-without][[l]]), "theta")
                                    thetas = .read.params(tree.nn); values = sort(.read.params(tree[-without][[l]]), decreasing = FALSE)
                                    expr = .constr.expr(tree.nn, type); leaves = .get.leaves(tree[-without][[l]]); d = length(leaves) + 1
                                    ff = .d.dell(parse(text=expr), c("leaf.new", leaves, thetas[1], thetas[-1][order(values)]), order = d)
                                    for(k in 1:(d-1)){formals(ff)[[k+1]]=X[,leaves[k]]}
                                    for(k in 1:length(values)){formals(ff)[[d+1+k]]=values[k]}
                                    ff.done = TRUE
                                }
                                coln = c(names[i], names[j])[which(c(names[i], names[j])!="tree")]
                                matr.p[i, j] = matr.p[j, i] = optimise(function(y){sum(log(attr(ff(leaf.new = X[, coln], theta = tau2theta(y, type)), "gradient")))}, interval = c(1e-8, theta2tau(select[3], type) - 1e-8), maximum = TRUE)$maximum
                        }}else{ 
                            if(is.null(without)){
                                tree.nn = list(.tree.without.params(c(tree, 1)))
                                values = sort(.read.params(c(tree, 1)), decreasing = FALSE)
                                leaves = .get.leaves(c(tree, 1))
                            }else{
                                tree.nn = list(.tree.without.params(c(tree[-without], 1)))
                                values = sort(.read.params(c(tree[-without], 1)), decreasing = FALSE)
                                leaves = .get.leaves(c(tree[-without], 1))
                            }
                           thetas = .read.params(tree.nn)
                           expr = .constr.expr(tree.nn, type)
                           d = length(leaves)
                           fp = .d.dell(parse(text=expr), c(leaves, thetas[order(values)]), order = d)
                           for(k in 1:d){formals(fp)[[k]]=X[,leaves[k]]}
                           for(k in 2:length(values)){formals(fp)[[d+k]]=values[k]}
                           matr.p[i, j] = matr.p[j, i] = optimise(function(y){sum(log(attr(fp(theta1.1 = tau2theta(y, type)), "gradient")))}, interval = c(1e-8, theta2tau(select[3], type) - 1e-8), maximum = TRUE)$maximum
                        }}
            }}
            
            select = c(min(row(matr.p)[which(matr.p==max(matr.p))]), max(col(matr.p)[which(matr.p==max(matr.p))]), tau2theta(max(matr.p), type))
            s = select[1:2]
            tree.n = list(c(tree[s], select[3]))
            
            if(class(tree.n[[1]][[length(tree.n[[1]])]])=="numeric"){tree.n = list(.union(tree.n[[1]], epsilon = epsilon, method = agg.method, ...))}
            
            tree = c(tree[-s], tree.n); names = c(names[-s], "tree")
            if(any(names!="tree")){without = which(names!="tree")}else{without = NULL}
            main.dim = main.dim - 1
        }
        .union(tree.n[[1]], epsilon = epsilon, method = agg.method, ...)
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
	if(is.null(margins) | (class(X)!="matrix") | (NROW(X)==1)){
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
		}else{
		.opt.margin(data = X, spec = spec)
		}
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
		op = optimise(f = function(par, data){sum(log(eval(do.call(paste("d", spec, sep = ""), args = list(x = data, par)))))}, data = data, lower = 0.0001, upper = 100, maximum = TRUE)$maximum
		eval(do.call(paste("p", spec, sep = ""), args = list(q = data, op)))	
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.max.min = function(X){
	if((any(X<=0)) |  (any(X>=1))){
		X[which(X >= 1)] = 1-1e-06
		X[which(X <= 0)] = 1e-06
		X}
	else{
		X}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.constraints.ui = function(tree, m, values){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
 
     if(any(s)){
         if(any(!s)){
            n.constr = length(which(!s))
            m.new = matrix(0, nrow = n.constr, ncol = length(values))
         	params = sapply(tree[which(!s)], function(r)r[[length(r)]])
			for(i in 1:n.constr){
         		m.new[i, which(values==params[i])]=1
         		m.new[i, which(values==tree[[n]])]=-1
         	}
            m = rbind(m, m.new)
            for(i in which(!s)){
            	m = .constraints.ui(tree[i], m, values)
            }
         }else{
            m = m
     }}else{
       m.new = matrix(0, nrow = (n-1), ncol = length(values))
       params = sapply(tree[-n], function(r)r[[length(r)]])
            for(i in 1:(n-1)){
         		m.new[i, which(values==params[i])]=1
         		m.new[i, which(values==tree[[n]])]=-1
         	}
        m = rbind(m, m.new)
            for(i in 1:(n-1)){
            	m = .constraints.ui(tree[i], m, values)
            }
    }
    return(m)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.rebuild = function(tree, values, theta){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     tree[[n]] = theta[which(values==tree[[n]])]
              
     if(any(s)){
         if(any(!s)){
            tree=c(tree[which(s)], lapply(tree[which(!s)], .rebuild, values, theta), tree[[n]])           
        }else{
            tree=tree
     }}else{
        tree = c(lapply(tree[-n], .rebuild, values, theta), tree[[n]])
     }
     return(tree)
}
