# estimate.r ############################################################################################################# 
# FUNCTION:               	DESCRIPTION: 
#  estimate.copula			    Estimates the structure and the parameter of a HAC for a given sample. 
#  .QML                     Estimation procedures based on binary trees and QML, i.e., for method = ML. (Internal function)
#  .PML                     Estimation procedures based on penalized QML, i.e., for method = PML. (Internal function)
#  .QML.hac                 Estimation procedures based on method = 1 and a prespecified hac-structure, i.e., hac != NULL. (Internal function)
#  .FML                     Full Maximum Likelihood (FML) estimation procedure. It needs an 'hac' object as argument to construct the log-likelihood which depends on the structure of the HAC. (Internal function)
#  .RML                     Recursive Maximum Likelihood (RML) estimation procedure. (Internal function)
#  .ub         			 	      Ensures the dependency parameter of the initial node being smaller than parameter of consecutive nodes. (Internal function) 
#  .margins				          Estimates the marginal distributions and returns the fitted values for a d-dimensional sample. (Internal function)   
#  .one.mar				          Estimates one marginal distributions for a given univariate sample. (Internal function)   
#  .max.min					        0's contained in the data matrix are set to 1e-16 and 1's to 1-1e-16. (Internal function) 
#  .constraints.ui          Returns a matrix of constraints according to the matrix ui of constrOptim. This matrix ensures the parameters being increasing from the highest to the lowest hierarchical level for the full ML approach. (Internal function)
#  .rebuild                 Matches the tree of a 'hac' object according to an ordered parameter vector. (Internal function) 
##########################################################################################################################

estimate.copula = function(X, type = 1, method = 1, hac = NULL, epsilon = 0, agg.method = "mean", margins = NULL, na.rm = FALSE, max.min = TRUE, ...){
	
	if(is.null(colnames(X))){g.names = names = paste0("X", 1:NCOL(X))}else{names = colnames(X)}
	
	X = .margins(X, margins)
	colnames(X) = names
	
	if(na.rm){X = na.omit(X, ...)}

	if(max.min){X = .max.min(X)}
	
	d = NCOL(X)	
    if(((type == 1) | (type == 3) | (type == 5) | (type == 7) | (type == 9)) & (d > 2)){
    	if(method == 1){
    	      if(is.null(hac)){
                res = .QML(X = X, type = type, epsilon = epsilon, agg.method = agg.method, names = names, ...)
            }else{
                res = .QML.fixed.tree(tree = hac$tree, X = X, type = hac$type)
        	  }    
        }else
        if(method == 2){
            if(is.null(hac)){
                stop("A hac object is required.")
            }else{
                res = .FML(X = X, type = type, hac = hac)
        	  }
        }else
        if(method == 3){
            res = .RML(X = X, type = type, method = method, epsilon = epsilon, agg.method = agg.method, names = names, ...)
        }else
        if(method == 4){
            if(type == 5 | type == 7){
              stop("Estimation of penalized HAC for this family is currently not supported.")
            }else{
              res = .PML(X = X, type = type, names = names)
            }
        }
    }else{
    	if((type == 2) | (type == 1)){
            res = c(as.list(names), fitCopula(gumbelCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
    	}else 
    	if((type == 4) | (type == 3)){
            res = c(as.list(names), fitCopula(claytonCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == 6) | (type == 5)){
            res = c(as.list(names), fitCopula(frankCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == 8) | (type == 7)){
            res = c(as.list(names), fitCopula(joeCopula(param = 1.5, dim = d), X, method = "ml")@estimate)
        }else 
    	if((type == 10) | (type == 9)){
            res = c(as.list(names), fitCopula(amhCopula(param = 0.25, dim = d), X, method = "ml")@estimate)
        }
    }
	hac(type = type, tree = res)
}     

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.QML = function(X, type, epsilon, agg.method = "mean", names, ...){
        main.dim = NCOL(X); tree = as.list(names);
        matr = matrix(0, main.dim, main.dim)
        
        upper.tau = if((type == 10) | (type == 9)){1/3 - 1e-10}else{1 - 1e-10}
        for(i in 1:(main.dim-1)){
            for(j in (i+1):main.dim){
                matr[i, j] = matr[j, i] = optimise(f = function(y, i, j){sum(log(.dAC(X[, i], X[, j], tau2theta(y, type), type)))}, i = i, j = j, interval = c(1e-08, upper.tau), maximum = TRUE)$maximum
            }
        }
        
        pair = c(min(row(matr)[which(matr == max(matr))]), max(col(matr)[which(matr == max(matr))]))
        current.theta = tau2theta(max(matr), type)

	    X[, pair[1]] = .cop.T(sample = X[, pair], theta = current.theta, type = type)
    		X = X[, -pair[2]]

	    colnames(X)[pair[1]] = "tree"
    		tree[[pair[1]]] = c(tree[pair], current.theta)
        tree = tree[-pair[2]]
		
		matr = matr[-pair[2],-pair[2]]
		Index.j = 1:(main.dim <- main.dim - 1); 
		
        while(main.dim >= 2){
        	pair = pair[1]
        	current.names = colnames(X)
            for(j in Index.j[-pair]){
                if ((current.names[pair] == "tree") & (current.names[j] != "tree")) {
                  upper.tau = theta2tau(tree[[pair]][[length(tree[[pair]])]], type)
                }
                else if ((current.names[pair] != "tree") & (current.names[j] == "tree")) {
                  upper.tau = theta2tau(tree[[j]][[length(tree[[j]])]], type)
                }
                else if ((current.names[pair] == "tree") & (current.names[j] == "tree")) {
                  upper.tau = min(theta2tau(c(tree[[pair]][[length(tree[[pair]])]], tree[[j]][[length(tree[[j]])]]), type))
                }
                matr[pair, j] = matr[j, pair] = optimise(f = function(y, j){sum(log(.dAC(X[, pair], X[, j], tau2theta(y, type), type)))}, j = j, interval = c(1e-08, upper.tau), maximum = TRUE)$maximum
            }
            
        		pair = c(min(row(matr)[which(matr == max(matr))]), max(col(matr)[which(matr == max(matr))]))
        		current.theta = tau2theta(max(matr), type)
        		
        		tree[[pair[1]]] = c(tree[pair], current.theta)
	        tree = tree[-pair[2]]

    	    		if((main.dim <- main.dim - 1) == 1){
    	    			return(.union(tree[[1]], epsilon = epsilon, method = agg.method, ...))
    	    		}
        		
		    X[, pair[1]] = .cop.T(sample = X[, pair], theta = current.theta, type = type)
    			X = X[, -pair[2]]
        	
		    colnames(X)[pair[1]] = "tree"

			matr = matr[-pair[2],-pair[2]]
			Index.j = 1:main.dim; 
		}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.PML = function(X, type, names){
  main.dim = NCOL(X); tree = as.list(names); n = NROW(X)
  matr.logL = matr = matrix(0, main.dim, main.dim)
  
  XX = X
  
  g.upper.tau = if((type == 10) | (type == 9)){
                    1/3 - 1e-08
                  }else{
                    1 - 1e-08}
  for(i in 1:(main.dim - 1)){
    for(j in (i + 1):main.dim){
      .min = optimise(f = function(y, i, j){-sum(log(.dAC(X[, i], X[, j], tau2theta(y, type), type)))}, i = i, j = j, interval = c(1e-08, g.upper.tau))
      matr[i, j] = matr[j, i] = .min$minimum
      matr.logL[i, j] = matr.logL[j, i] = .min$objective
    }
  }
  
  pair = c(min(row(matr.logL)[which(matr.logL == min(matr.logL))]), max(col(matr.logL)[which(matr.logL == min(matr.logL))]))
  current.theta = tau2theta(matr[pair[1], pair[2]], type)
  
  X[, pair[1]] = .cop.T(sample = X[, pair], theta = current.theta, type = type)
  X = X[, -pair[2]]
  
  tree[[pair[1]]] = c(tree[pair], current.theta)
  tree = tree[-pair[2]]
  
  colnames(X)[pair[1]] = paste(.get.leaves(tree[[pair[1]]]), collapse = "")
  
  matr = matr[-pair[2],-pair[2]]; matr.logL = matr.logL[-pair[2],-pair[2]]
  Index.j = 1:(main.dim <- main.dim - 1)
  
  while(main.dim >= 2){
    pair = pair[1]
    current.names = colnames(X)
    for(j in Index.j[-pair]){
      if(!(current.names[pair] %in% names) & (current.names[j] %in% names)) {
        upper.tau = theta2tau(tree[[pair]][[length(tree[[pair]])]], type)
      }
      else if ((current.names[pair] %in% names) & !(current.names[pair] %in% names)) {
        upper.tau = theta2tau(tree[[j]][[length(tree[[j]])]], type)
      }
      else if (!(current.names[pair] %in% names) & !(current.names[pair] %in% names)) {
        upper.tau = min(theta2tau(c(tree[[pair]][[length(tree[[pair]])]], tree[[j]][[length(tree[[j]])]]), type))
      }
      .min = optimise(f = function(y, j){-sum(log(.dAC(X[, pair], X[, j], tau2theta(y, type), type)))}, j = j, interval = c(1e-08, upper.tau))
      matr[pair, j] = matr[j, pair] = .min$minimum
      matr.logL[j, pair] = matr.logL[pair, j] =.min$objective   
    }
    
    pair = c(min(row(matr.logL)[which(matr.logL == min(matr.logL))]), max(col(matr.logL)[which(matr.logL == min(matr.logL))]))
    current.theta = tau2theta(matr[pair[1], pair[2]], type)
    
    sub.tree = c(tree[pair], current.theta)
    sub.d = length(sub.tree)
    s = sapply(sub.tree, is.character)[-sub.d]
    sub.X = X[, c(unlist(sub.tree[which(s)]), unlist(sapply(sub.tree[which(!s)], FUN = function(r){paste(.get.leaves(r), collapse = "")})))]
    
    repeat{
      n.pars = length(.read.params(sub.tree))
      if(n.pars > 1){
        .trees = which(sapply(sub.tree, is.list))
        
        bic.full = numeric(length(.trees))
        for(k in 1:length(.trees)){
          .sub.tree = sub.tree
          .sub.sub.tree = .sub.tree[[.trees[k]]]
          .sub.d = length(.sub.tree)
          
          .s = sapply(.sub.tree, is.character)[-sub.d]
          
          .s.k = sapply(.sub.sub.tree, is.character)[-length(.sub.sub.tree)]
          .sub.names = unlist(.sub.sub.tree[which(.s.k)])
          if(any(.s)){
            .r.names = unlist(.sub.tree[which(.s)])
            .X = cbind(XX[, c(.r.names, .sub.names)])
            if(is.null(.sub.names)){colnames(.X) = .r.names}
          }else{
            .X = cbind(XX[, .sub.names])
            if(length(.sub.names) == 1){colnames(.X) = .sub.names}
          }
          if(any(!.s.k)){
            for(nk in which(!.s.k)){
              pseudo.var = .cop.transform(XX[, .get.leaves(.sub.sub.tree[[nk]])], .sub.sub.tree[[nk]], type = type)
              .sub.tree[[.trees[k]]][[nk]] = paste(.get.leaves(.sub.sub.tree[[nk]]), collapse = "")
              .X = cbind(.X, pseudo.var)
              colnames(.X)[NCOL(.X)] = .sub.tree[[.trees[k]]][[nk]]
            }
          }
          
          if(length(.trees) > 1){
            for(l in (1:length(.trees))[-k]){
              pseudo.var = .cop.transform(XX[, .get.leaves(.sub.tree[[.trees[l]]])], .sub.tree[[.trees[l]]], type = type)
              .sub.tree[[.trees[l]]] = paste(.get.leaves(.sub.tree[[.trees[l]]]), collapse = "")
              .X = cbind(.X, pseudo.var)
              colnames(.X)[NCOL(.X)] = .sub.tree[[.trees[l]]] 
            }
          }
          bic.full[k] = n.pars * log(n) - 2 * .logLik_oneHierachy(U = .X, tree = .sub.tree, type = type)  
        }
        
        bic.restr = numeric(length(.trees))
        for(k in 1:length(.trees)){
          .sub.tree = sub.tree
          .sub.sub.tree = .sub.tree[[.trees[k]]]
          .sub.d = length(.sub.tree)
          
          .s = sapply(.sub.tree, is.character)[-sub.d]
          
          .s.k = sapply(.sub.sub.tree, is.character)[-length(.sub.sub.tree)]
          .sub.names = unlist(.sub.sub.tree[which(.s.k)])
          if(any(.s)){
            .X = XX[, c(unlist(.sub.tree[which(.s)]), .sub.names)]
          }else{
            .X = XX[, .sub.names]
          }
          if(any(!.s.k)){
            for(nk in which(!.s.k)){
              pseudo.var = .cop.transform(XX[, .get.leaves(.sub.sub.tree[[nk]])], .sub.sub.tree[[nk]], type = type)
              .X = cbind(.X, pseudo.var)
            }
          }
          
          upper.tau = g.upper.tau
          if(length(.trees) > 1){
            for(l in (1:length(.trees))[-k]){
              upper.tau = c(upper.tau, theta2tau(.sub.tree[[.trees[l]]][[length(.sub.tree[[.trees[l]]])]], type = type))
              pseudo.var = .cop.transform(XX[, .get.leaves(.sub.tree[[.trees[l]]])], .sub.tree[[.trees[l]]], type = type)
              .X = cbind(.X, pseudo.var)
            }
          }
          upper.tau = min(upper.tau)
          bic.restr[k] = optimise(f = function(tau){
                                        (n.pars - 1) * log(n) - 2 * sum(.d.multi.AC(.X, tau2theta(tau, type), type + 1))
                                     }, interval = c(0, upper.tau))$objective
        }
        
        if(min(bic.restr - bic.full) < 0){
          .tree = .trees[which.min(bic.restr - bic.full)]
          sub.sub.tree = sub.tree[[.tree]]
          sub.tree = c(sub.sub.tree[-length(sub.sub.tree)], sub.tree[-.tree])
          
          sub.d = length(sub.tree)
          sub.tree[[sub.d]] = Inf
          
          s = sapply(sub.tree, is.character)[-sub.d]
          sub.names = sapply(sub.tree[which(!s)], FUN = function(r){paste(.get.leaves(r), collapse = "")})
          if(length(sub.names) == 0){
            sub.X = XX[, unlist(sub.tree[which(s)])]
          }else{
            sub.X = XX[, unlist(sub.tree[which(s)])]
            for(l in which(!s)){
              sub.X = cbind(sub.X, .cop.transform(XX[,.get.leaves(sub.tree[[l]])], sub.tree[[l]], type))
            }
            colnames(sub.X) = c(unlist(sub.tree[which(s)]), sub.names)
          }
          
          if(any(!s)){
            upper.tau = theta2tau(min(unlist(sapply(sub.tree[which(!s)], function(r){r[length(r)]}))), type = type)  
          }else{
            upper.tau = g.upper.tau
          }
          
          current.theta = sub.tree[[sub.d]] = tau2theta(optimise(f = function(y){-sum(.d.multi.AC(sub.X, tau2theta(y, type), type + 1))}, interval = c(1e-08, upper.tau))$minimum, type)
        }else{break}
      }else{break}
    }
    
    tree[[pair[1]]] = sub.tree
    tree = tree[-pair[2]]
    
    main.dim <- main.dim - 1
    if(main.dim == 1){
        return(tree[[1]])
    }
    
    X[, pair[1]] = .cop.transform(XX[,.get.leaves(sub.tree)], sub.tree, type)
    X = X[, -pair[2]]
    colnames(X)[pair[1]] = paste(.get.leaves(sub.tree), collapse = "")
    
    matr = matr[-pair[2], -pair[2]]; matr.logL = matr.logL[-pair[2], -pair[2]];
    Index.j = 1:main.dim
  }
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.QML.fixed.tree = function(tree, X, type){
          if(length(tree)==1){tree = tree[[1]]}
	      n = length(tree);
	      s = sapply(tree, is.character)

	      if(any(!s[-n])){
                Tau = NULL
			    for(j in which(!s[-n])){
                    .names.j = .get.leaves(tree[[j]])
                    tree[[j]] = .QML.fixed.tree(tree[[j]], X[,.names.j], type)
                    X = cbind(X, .cop.transform(X[,.names.j], tree[[j]], type))
                    X = X[,-which(colnames(X) %in% .names.j)]; colnames(X) = c(colnames(X)[-NCOL(X)], paste("tree", j, sep = ""))
                    Tau = c(Tau, theta2tau(tree[[j]][[length(tree[[j]])]], type = type))
                }
                tree[[n]] = tau2theta(optimise(f = function(y){sum(.d.multi.AC(X, tau2theta(y, type), type + 1))}, interval = c(1e-08, min(Tau)), maximum = TRUE)$maximum, type)
		    }else{
		        tree[[n]] = estimate.copula(X, type = type + 1)$tree[[NCOL(X)+1]]
	      }
	      tree
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.FML = function(X, type, hac){
    values = get.params(hac, sort.v = TRUE, decreasing=FALSE)
    tree.full = hac$tree
	initial=if((type == 1) | (type == 7)){1+1e-8}else{1e-8}
    ui = .constraints.ui(tree.full, m = matrix(c(1, rep(0, length(values)-1)), nrow=1), values = values)
    LL = to.logLik(X, hac)
    optim = constrOptim(values, f=LL, grad=NULL, ui=as.matrix(ui), ci=as.vector(c(initial, rep(1e-8, NROW(ui)-1))), control=list(fnscale=-1), hessian=FALSE)
    .rebuild(tree.full, values, optim$par)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.RML = function(X, type, method, epsilon, agg.method, names, ...){
    main.dim = NCOL(X); tree = as.list(names)
    select = if((type == 10) | (type == 9)){c(0, 0, 1-1e-8)}else{c(0, 0, 100)}
        
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
            tree.n = c(tree[s], select[3])
            
            if(is.numeric(tree.n[[length(tree.n)]])){tree.n = .union(tree.n, epsilon = epsilon, method = agg.method, ...); select[3] = tree.n[[length(tree.n)]]}
            
            tree = c(tree[-s], list(tree.n)); names = c(names[-s], "tree")
            if(any(names!="tree")){without = which(names!="tree")}else{without = NULL}
            main.dim = main.dim - 1
        }
        tree.n
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.ub = function(tree.1, tree.2, type){
  	if((is.numeric(tree.1)) & is.character(tree.2))
  		theta2tau(tree.1, type)
  	else 
  	if(is.character(tree.1) & is.character(tree.2))
  		theta2tau(tree.2, type)
	  else 
	  if(is.numeric(tree.1) & is.numeric(tree.2))
  		theta2tau(min(tree.1, tree.2), type)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.margins = function(X, margins, ...){
	if(is.null(margins) | !inherits(X, "matrix") | (NROW(X) == 1)){
		X
	}else{
		if(length(margins)==1){
		X = apply(X, 2, .one.mar, spec = margins, ...)
	}else{
		for(i in 1:NCOL(X)){X[,i] = .one.mar(X[,i], margins[i], ...)}}
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
		op = constrOptim(theta = c(1, 1), f = loglik, grad = NULL, ui = matrix(c(1, 0, -1, 0, 0, 1), nrow = 3, byrow = TRUE), ci = c(-rep(boundary, 2), 0), data = data, control = list(fnscale = -1), hessian = FALSE)
		eval(do.call(paste("p", spec, sep = ""), args = list(q = data, op$par[1], op$par[2])))					
	}else{
	if((spec == "exp")){
		op = optimise(f = function(par, data){sum(log(eval(do.call(paste("d", spec, sep = ""), args = list(x = data, par)))))}, data = data, lower = 0.0001, upper = 100, maximum = TRUE)$maximum
		eval(do.call(paste("p", spec, sep = ""), args = list(q = data, op)))	
	}}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.max.min = function(X){
	if((any(X <= 0)) |  (any(X >= 1))){
		X[which(X >= 1)] = 1-1e-16
		X[which(X <= 0)] = 1e-16
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
