# copula_pdf.r ###########################################################################################################
# FUNCTION:               	DESCRIPTION:
#  dAC							Computes the values of the density of 2-dimensional copulae.
#  .gumb.12.density		2-dim density of Gumbel copulae. (Internal function)
#  .clay.12.density		2-dim density of Clayton copulae. (Internal function)
#  dHAC						Computes the values of the density of 2- or 3-dimensional copulae.
#  .gumb.12.3.density   3-dim density of HAC Gumbel copulae. (Internal function)
#  .clay.12.3.density  	 3-dim density of HAC Clayton copulae. (Internal function)     
##########################################################################################################################

dAC = function(x, y, theta = 1.0, type = AC_GUMBEL, na.rm = FALSE, max.min = TRUE){
		
	if(max.min == TRUE){
		x = .max.min(x)
		y = .max.min(y)}
	
	if(na.rm == TRUE){
		x = .na.rm(x)
		y = .na.rm(y)}
	
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

dHAC = function(X, hac, margins = NULL, na.rm = FALSE, max.min = TRUE){
	X = .margins(X, margins)
	
	if(max.min == TRUE){
			X = .max.min(X)}
			
	if(na.rm == TRUE){
			X = .na.rm(X)}
	
	cop.12.3.density = function(u1, u2, u3, theta1, theta2, type){
		if(type == HAC_GUMBEL)
			.gumb.12.3.density(u1, u2, u3, theta1, theta2)
		else if(type == HAC_CLAYTON)
			.clay.12.3.density(u1, u2, u3, theta1, theta2)
		else if(type == HAC_ROTATED_GUMBEL)
			.gumb.12.3.density(1 - u1, 1 - u2, 1 - u3, theta1, theta2)
	}
	
    if((hac$type == HAC_GUMBEL) || (hac$type == HAC_CLAYTON) || (hac$type == HAC_ROTATED_GUMBEL)){
        if(dim(X)[2] != 3){
            return("The PDF of HAC for d>3 is not implemented yet")}
    else{
    if(class(hac$model$V1) == "character"){
        a = X[,hac$model$V1]
        X = X[, -which(hac$model$V1 == colnames(X))]
        X = cbind(X, a)
        l.theta = c(hac$model$theta, hac$model$V2$theta)}
    else 
    if(class(hac$model$V2) == "character"){
        a = X[,hac$model$V2]
        X = X[, -which(hac$model$V2 == colnames(X))]
        X = cbind(X, a)
        l.theta = c(hac$model$theta, hac$model$V1$theta)}
    return(cop.12.3.density(X[,1], X[,2], X[,3], l.theta[1], l.theta[2], hac$type))}
    }else if(hac$type == GAUSS){
        return(dcopula(normalCopula(hac$model[lower.tri(hac$model)], dim = dim(X)[2], dispstr = "un"), X))
    }else if(hac$type == AC_GUMBEL){
        return(dcopula(gumbelCopula(hac$model$theta, dim = hac$model$dim), X))
    }else if(hac$type == AC_CLAYTON){
        return(dcopula(claytonCopula(hac$model$theta, dim = hac$model$dim), X))
    }
}

#-------------------------------------------------------------------------------------------------------------------------------

	.gumb.12.3.density = function(u1, u2, u3, theta1, theta2){
##  theta1 < theta2 ##
		l1 = -log(u1)
		l2 = -log(u2)
		l3 = -log(u3)
		theta1m1 = (theta1 - 1)
		theta2m1 = (theta2 - 1)
		onedtheta1 = 1/theta1
		l12.2 = (l1)^theta2 + (l2)^theta2
		l12.2.t2 = l12.2^(1/theta2)
		c12 = (l12.2.t2)^theta1
		c12.3 = c12 + l3^theta1
		tc12.3 = (theta1m1 + c12.3^onedtheta1)
		
		-((l12.2.t2^(theta1 - 2)*(l1^theta2m1)*l12.2^(-2 + 1/theta2)*l2^theta2m1*(theta2m1*(-l12.2.t2) * tc12.3 * c12.3 + l12.2.t2*(-(c12 * (theta1m1 * theta1 + 2 * theta1m1 * c12.3^onedtheta1 + c12.3^(2*onedtheta1))) + theta1m1*tc12.3 * l3^theta1))*c12.3^(onedtheta1-3)*l3^theta1m1) / (exp(c12.3^onedtheta1)*u1*u2*u3))
	}
	
#-------------------------------------------------------------------------------------------------------------------------------
	
	.clay.12.3.density = function(u1,u2,u3,theta2,theta1){
		u1t1u2t1 = -1+u1^(-theta1)+u2^(-theta1)
		it1 = 1/theta1
		it2 = 1/theta2
		u1mt1 = u1^(-1-theta1)
		u2mt1 = u2^(-1-theta1)
		u3mt2 = u3^(-1-theta2)
		u1t1u2t1it1 = u1t1u2t1^(-it1)
		shorter = -1+u1t1u2t1it1^(-theta2)+u3^(-theta2)
		long = shorter^(-2-it2)
		
		(-2-it2)*(-1-it2)*(theta2^2)*u1mt1*u2mt1*(u1t1u2t1^(-2-2*it1))*(u1t1u2t1it1^(-2-2*theta2))*u3mt2*(shorter^(-3-it2))-(-1-it2)*(-1-theta2)*theta2*u1mt1*u2mt1*(u1t1u2t1^(-2-2/theta1))*(u1t1u2t1it1^(-2-theta2))*u3mt2*long+(-1-it1)*theta1*(-1-it2)*theta2*u1mt1*u2mt1*(u1t1u2t1^(-2-it1))*(u1t1u2t1it1^(-1-theta2))*u3mt2*long
	}