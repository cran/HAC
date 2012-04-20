# copula_pdf.r ###########################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .dAC						Computes the values of the density of 2-dimensional copulae. (Internal function)
#  .gumb.12.density			2-dim density of Gumbel copulae. (Internal function)
#  .clay.12.density			2-dim density of Clayton copulae. (Internal function)
#  dHAC						Computes the values of the density of 2- or 3-dimensional copulae.
#  .cop.pdf					Computes the values of the density. (Internal function)
#  .d.dell                  Computes the deriavtive of an expression given by .constr.expr. (Internal function)
#  .constr.expr             A recursive function. Computes an expression of the cdf for a given type of copula.(Internal function)
##########################################################################################################################

.dAC = function(x, y, theta = 1.0, type = AC_GUMBEL){	
	if((type == HAC_GUMBEL) || (type == AC_GUMBEL)){
		.gumb.12.density(x, y, theta)
	}else if((type == HAC_CLAYTON) || (type == AC_CLAYTON)){
		.clay.12.density(x, y, theta)       
	}else if(type == GAUSS){
		dcopula(normalCopula(theta, 2, dispstr = "un"), cbind(x, y))
	}else if(type == HAC_ROTATED_GUMBEL){
		.gumb.12.density(1 - x, 1 - y, theta)
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
	X = .margins(X, margins)
			
	if(na.rm == TRUE){
			X = na.omit(X, ...)}
		
    if((hac$type == HAC_GUMBEL) | (hac$type == HAC_CLAYTON)){
        if(NCOL(X) >= 3){
        	.cop.pdf(tree = hac$tree, sample = X, type = hac$type, eval = eval)
    }}else if(hac$type == GAUSS){
        return(dcopula(normalCopula(hac$tree[lower.tri(hac$model)], dim = dim(X)[2], dispstr = "un"), X))
    }else if(hac$type == AC_GUMBEL){
    	n = length(hac$tree)
        return(dcopula(gumbelCopula(hac$tree[[n]], dim = (n-1)), X))
    }else if(hac$type == AC_CLAYTON){
    	n = length(hac$tree)
        return(dcopula(claytonCopula(hac$tree[[n]], dim = (n-1)), X))
    }else if(hac$type == HAC_ROTATED_GUMBEL){
        stop("HAC-pdf for HAC_ROTATED_GUMBEL is not implemented yet.")
    }
}

#-------------------------------------------------------------------------------------------------------------------------------

.cop.pdf = function(tree, sample, type, eval){
	d = NCOL(sample); names = colnames(sample)
	string.expr = .constr.expr(tree, type)
    Dd = .d.dell(parse(text=string.expr), names, d)
    
    if(eval){
        for(i in 1:d){formals(Dd)[[i]]=sample[,i]}
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
 
     if(any(s==TRUE)){
         if(any(s==FALSE)){
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
