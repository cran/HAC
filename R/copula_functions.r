# copula_functions.r #####################################################################################################
# FUNCTION:               DESCRIPTION:
#  definitions					Global paramters.
#  theta.eps					Is addressed during the estimation procedure.
#  TAU							Possible argument of estimate.copula.
#  ML								Possible argument of estimate.copula.
#  tau2theta       				Convertes Kendall's rank correlation coefficient into dependence parameter.
#  theta2tau          			Convertes dependence parameter into Kendall's rank correlation coefficient.
#  phi        						The generator function.
#  phi.inv				  		The inverse of the generator function.
#  copMult        				Computes the value of d-dimensional AC.      
##########################################################################################################################

HAC_GUMBEL = 0
AC_GUMBEL = 1
HAC_ROTATED_GUMBEL = 2
HAC_CLAYTON = 3
AC_CLAYTON = 4
GAUSS = 5

TAU = 0
ML  = 1

#-------------------------------------------------------------------------------------------------------------------------------

tau2theta = function(tau, type = HAC_GUMBEL){
	n = length(tau)
	if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == HAC_ROTATED_GUMBEL) || (type == AC_CLAYTON) || (type == HAC_CLAYTON)){
		for(i in 1 : n){if((tau[i] < 0) | (tau[i] > 1)){return(warning(paste("Element[", i,"] should be in [0, 1).")))}}
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == HAC_ROTATED_GUMBEL))
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
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == HAC_ROTATED_GUMBEL)){
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
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == HAC_ROTATED_GUMBEL)){
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
	
    if((type == AC_GUMBEL) || (type == HAC_GUMBEL) || (type == HAC_ROTATED_GUMBEL)){
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
