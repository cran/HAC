# copula_pdf.r ###########################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .dAC						Computes the values of the density of 2-dimensional copulae. (Internal function)
#  .gumb.12.density			2-dim density of Gumbel copulae. (Internal function)
#  .clay.12.density			2-dim density of Clayton copulae. (Internal function)
#  dHAC						Computes the values of the density.
#  .cop.pdf					Supplementary function of dHAC. (Internal function)
#  .d.dell                  Computes the deriavtive of an expression given by .constr.expr. (Internal function)
#  .constr.expr             Computes an expression of the cdf for a given copula type. (Internal function)
#  to.logLik                Returns the log-Likelihood function or its value.
#  .tree.without.params     Tranforms a tree of a hac object to a tree with symbolic parameter. (Internal function)
##########################################################################################################################

.dAC = function(x, y, theta = 1.0, type = AC_GUMBEL){	
	if((type == HAC_GUMBEL) | (type == AC_GUMBEL)){
		.gumb.12.density(x, y, theta)
	}else if((type == HAC_CLAYTON) | (type == AC_CLAYTON)){
		.clay.12.density(x, y, theta)       
	}else if(type == GAUSS){
		dcopula(normalCopula(theta, 2, dispstr = "un"), cbind(x, y))
	}
}

#-------------------------------------------------------------------------------------------------------------------------------

.gumb.12.density = function(x, y, theta){
	lu1 = -log(x)
	lu2 = -log(y)
	(lu1^(-1 + theta)*(-1 + theta + (lu1^theta + lu2^theta)^(1/theta))*(lu1^theta + lu2^theta)^(-2 + 1/theta)*lu2^(-1 + theta))/(exp((lu1^theta + lu2^theta)^(1/theta))*x*y)
}
	
#-------------------------------------------------------------------------------------------------------------------------------
	
.clay.12.density = function(x, y, theta){
	u1pt = x^(-theta)
	u2pt = y^(-theta)
	(1+theta)*u1pt*u2pt*((u1pt+u2pt-1)^(-1/theta - 2))/(x*y)
}
	
#-------------------------------------------------------------------------------------------------------------------------------

dHAC = function(X, hac, eval = TRUE, margins = NULL, na.rm = FALSE, ...){ 

    X = .one.ob(X, margins)
    names = colnames(X)    
    if(any(!(names %in% .get.leaves(hac$tree)))){stop("The colnames of X have to coincide with the specifications of the copula model hac.")}
			
	if(na.rm){X = na.omit(X, ...)}
    
    type = hac$type; d = NCOL(X)
    if(type != GAUSS){
        if((d >= 3) & ((type == HAC_GUMBEL) | (type == HAC_CLAYTON))){
           return(.cop.pdf(tree = hac$tree, sample = X, type = type, d = d, names = names, eval = eval))
        }else{
            colnames(X) = c()
            if(type == AC_GUMBEL){
                return(dcopula(gumbelCopula(hac$tree[[d+1]], dim = d), X)[-1])
            }else{
                if(type == AC_CLAYTON){
                    return(dcopula(claytonCopula(hac$tree[[d+1]], dim = d), X)[-1])
            }}
    }}else{ 
        colnames(X) = c()
        return(dcopula(normalCopula(hac$tree[lower.tri(hac$model)], dim = d, dispstr = "un"), X)[-1]) 
    }
}

#-------------------------------------------------------------------------------------------------------------------------------

.cop.pdf = function(tree, sample, type, d, names, eval){     
	string.expr = .constr.expr(tree, type)
    Dd = .d.dell(parse(text=string.expr), names, d)
    
    if(eval){
        for(i in 1:d){formals(Dd)[[i]]=sample[-1 ,i]}
        c(attr(Dd(), "gradient"))
    }else{
        Dd
    }
}

#-------------------------------------------------------------------------------------------------------------------------------

.d.dell = function(expr, name, order){
   if(order==1){
        deriv(expr, name[order], function.arg = name)
   }else{
        .d.dell(D(expr, name[order]), name, order-1)}
}

#---------------------------------------------------------------------------------------------------

.constr.expr = function(tree, type){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
 
     if(any(s)){
         if(any(!s)){
           if(type==HAC_GUMBEL){
                 paste("exp(-(", paste("(-log(", unlist(tree[which(s)]),"))^", tree[[n]], collapse="+", sep = ""),"+", paste("(-log(",sapply(tree[which(!s)], .constr.expr, type=type),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")
             }else{
                 paste("(", paste("(", unlist(tree[which(s)]),"^(-", tree[[n]],")-1)", collapse="+", sep = ""),"+", paste("((", sapply(tree[which(!s)], .constr.expr, type=type),")^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }
 }else{
             if(type==HAC_GUMBEL){
                 paste("exp(-(", paste("(-log(", unlist(tree[-n]),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")
             }else{
                 paste("(", paste("(",unlist(tree[-n]),"^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }
 }}else{
             if(type==HAC_GUMBEL){
                 paste("exp(-(", paste("(-log(", sapply(tree[-n], .constr.expr, type=type),"))^", tree[[n]], collapse="+", sep = ""),")^(1/", tree[[n]],"))", sep="")           
             }else{
                 paste("(", paste("((", sapply(tree[-n], .constr.expr, type=type),")^(-", tree[[n]],")-1)", collapse="+", sep = ""), "+1)^(-1/", tree[[n]], ")", sep="")
             }
}}

#---------------------------------------------------------------------------------------------------

to.logLik = function(X, hac, eval = FALSE, margins = NULL, na.rm = FALSE, ...){
	X = .margins(X, margins)
			
	if(na.rm){X = na.omit(X, ...)}
    
    tree = .tree.without.params(hac$tree)
    thetas = .read.params(tree); values = get.params(hac); d = NCOL(X)
    expr = .constr.expr(tree, hac$type)
    f = .d.dell(parse(text=expr), c(colnames(X), thetas[order(values)]), order=d)
    for(i in 1:d){formals(f)[[i]]=X[,i]}
    
    g = function(theta, density=f){
            n.par = length(theta)
            for(i in 1:n.par){formals(density)[[length(formals(density))-n.par+i]]=theta[i]}
            sum(log(c(attr(density(), "gradient"))))    
    }
        
    if(!eval){g}else{g(values)}
}
 
#---------------------------------------------------------------------------------------------------
 
.tree.without.params = function(tree, k=1, l=1){
     if(length(tree)==1){tree=tree[[1]]}
     n = length(tree)
     s = sapply(tree[-n], is.character)
     tree[[n]] = paste("theta",k,".",l, sep="")
     
     if(any(s)){
         if(any(!s)){
            for(i in which(!s)){
                tree[[i]]=.tree.without.params(tree[[i]], k=k+1,l=i)
            }       
         }else{
            tree = tree
         }}else{
         for(i in 1:(n-1)){
                tree[[i]]=.tree.without.params(tree[[i]], k=k+1,l=i)
         }}
    return(tree)        
}
