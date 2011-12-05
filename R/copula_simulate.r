# copula_simulate.r ######################################################################################################
# FUNCTION:               	DESCRIPTION:
#  rAC						Simulates values of AC.
#  .f_gumbel				Samples from the inverse Laplace-Stietjes transfrom of a Gumbel copula. (Internal function)
#  .f_clayton				Samples from the inverse Laplace-Stietjes transfrom of a Clayton copula. (Internal function)
#  rHAC					Simulates values of HAC.
#  .theta      				Computes the ratio of two dependency parameters. (Internal function)
#  .initial  					Samples from the inverse Laplace-Stietjes transfrom of the initial node of nested AC. (Internal function) 
#  .stayStage  			Simulates values of the initial node of HAC. (Internal function)     
#  .fReject 				The fast rejection algorithm is implemented. Have a look at Hofert (2010) "Efficiently sampling nested ...". (Internal function)     
#  .follow					Samples the inverse Laplace-Stietjes transfrom of all successive nodes of nested AC. (Internal function)   
#  .nextStage  			Simulates values of successive nodes of HAC. (Internal function)     
#  .simualte.hac  		Simulates all successive nodes. (Internal function)     
#  .rHAC 					Simulates the initial node. (Internal function)     
##########################################################################################################################

rAC = function(n, theta = 1.5, dim = 2, type = AC_GUMBEL){
    if(type == AC_GUMBEL)
	I = .f_gumbel(n, theta)
	else if(type == AC_CLAYTON)
	I = .f_clayton(n, theta)
	
    X = matrix(runif(dim * n), n, dim)
    phi(-log(X) / I, theta, type)
}

#-------------------------------------------------------------------------------------------------------------------------------

.f_gumbel = function(n, theta){
	as.vector(stabledist::rstable(n, 1 / theta, 1, cos(pi / (2 * theta))^(theta), 0, pm = 1))
}

#-------------------------------------------------------------------------------------------------------------------------------

.f_clayton = function(n, theta){
	as.vector(rgamma(n, shape = 1 / theta))
}

#-------------------------------------------------------------------------------------------------------------------------------

rHAC = function(n, hac){
	if(hac$type == AC_CLAYTON){
		res = rAC(n, theta = hac$model$theta, dim = hac$model$dim, type = AC_CLAYTON)
	}
	else
	if(hac$type == AC_GUMBEL){
		res = rAC(n, theta = hac$model$theta, dim = hac$model$dim, type = AC_GUMBEL)
	}
	else
	if(hac$type != HAC_ROTATED_GUMBEL){
		res = .rHAC(n, hac)
	}
	else
	if(hac$type == HAC_ROTATED_GUMBEL){
		hac$type = HAC_GUMBEL
		res = (1 - .rHAC(n, hac))
	}
	res
}

#-------------------------------------------------------------------------------------------------------------------------------

.theta = function(theta.i, theta.j){theta.i/theta.j}

#-------------------------------------------------------------------------------------------------------------------------------

.initial = function(n, Ltheta, type){
	mat = matrix(0, nrow = n)
		if(type == HAC_GUMBEL){
			mat = matrix(stabledist::rstable(n, alpha = 1 / Ltheta, beta = 1, gamma = cos(pi/(2 * Ltheta))^(Ltheta), delta = (Ltheta == 1) * 1, pm = 1), nrow = n)
			mat}
		else{
		if(type == HAC_CLAYTON){
			mat = matrix((rgamma(n, 1/Ltheta)), nrow = n)
			mat}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.fReject = function(alpha, I){
	if(alpha == 1){I}
	else{
		m = round(I) + (round(I) == 0) * 1
		m.max = max(m)
		m.length = length(m)
		
		gamma = (cos(alpha * pi / 2) * I * 1 / m)^(1 / alpha)	
		M = matrix(0, nrow = m.length, ncol = (m.max))
		
		for(i in 1:m.length){M[i, (1:m[i])] = 1}	
		M[which(M > 0)] = stabledist::rstable(n = sum(m), alpha = alpha, beta = 1, gamma = 1, delta = 0, pm = 1)
		G = matrix(rep(gamma, m.max), ncol = m.max)
		N =  M * G
		U = runif(m.length * m.max)
		Com = (U <= exp(-N))
		while(any((Com == FALSE) | (is.na(Com)))){
			fa = which((Com == FALSE) | (is.na(Com)))
			n = length(fa)
			N[fa] = stabledist::rstable(n, alpha = alpha, beta = 1, pm = 1) * G[fa]
			Com[fa] = (runif(n) <= exp(-N[fa]))}
		rowSums(N)}
}
												
#-------------------------------------------------------------------------------------------------------------------------------
										
.follow = function(n, Ltheta, I, type){
	mat = matrix(0, nrow = n)
		if(type == HAC_GUMBEL){
			gamma = (cos(pi/(2 * Ltheta)) * I)^(Ltheta)
			delta = (Ltheta == 1) * I
			mat = matrix(stabledist::rstable(n, alpha = 1 / Ltheta, beta = 1, pm = 1) * gamma + delta, nrow = n)
			mat}
		else{
		if(type == HAC_CLAYTON){
			mat = matrix(.fReject(alpha = 1 / Ltheta, I = I), nrow = n)
			mat}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.nextStage = function(n, d, Y, Ltheta, type){
	LU = matrix(rexp(n * d, rate = 1), nrow = n, ncol = d)
		if(type == HAC_GUMBEL){
			L =	phi(LU / matrix(c(rep(Y, d)), nrow = n), Ltheta, type = HAC_GUMBEL)
			L}
		else{
		if(type == HAC_CLAYTON){
			L =	phi(LU / matrix(c(rep(Y, d)), nrow = n), Ltheta, type = HAC_CLAYTON)
			L}}
}

#-------------------------------------------------------------------------------------------------------------------------------

.stayStage = function(n, d, Y, Ltheta, type){
	LU = matrix(rexp(n * d, rate = 1), nrow = n, ncol = d)
		if((type == HAC_GUMBEL) || (type == AC_GUMBEL)){
			L = phi(LU / matrix(c(rep(Y, d)), nrow = n), Ltheta, type = HAC_GUMBEL)
			L}
		else{
		if((type == HAC_CLAYTON) || (type == AC_CLAYTON)){
			L  = phi(LU / matrix(c(rep(Y, d)), nrow = n), Ltheta, type = HAC_CLAYTON)
			L}}
}
	
#-------------------------------------------------------------------------------------------------------------------------------

.simulate.hac = function(n, Lhac, First, ober.theta, type){
   	if((class(Lhac$V1) == "character") & (class(Lhac$V2) == "character")){
		Y1 = .follow(n, Ltheta = .theta(Lhac$theta, ober.theta), I = First, type)
		dim = length(Lhac$V1) + length(Lhac$V2)
		v = .nextStage(n, d = dim, Y = Y1, Ltheta = Lhac$theta, type)
		colnames(v) = c(Lhac$V1, Lhac$V2)}
    else if((class(Lhac$V1) != "character") & (class(Lhac$V2) == "character")){
    	dim = length(Lhac$V2)
    	Y2 = .follow(n, Ltheta = .theta(Lhac$theta, ober.theta), I = First, type = type)
        v1 = .simulate.hac(n, Lhac = Lhac$V1, First = Y2, ober.theta = Lhac$theta, type = type)
        v2 = .nextStage(n, d = dim, Y = Y2, Ltheta = Lhac$theta, type = type)
        colnames(v2) = Lhac$V2
        v = cbind(v1, v2)}
    else if((class(Lhac$V1) == "character") & (class(Lhac$V2) != "character")){
    	dim = length(Lhac$V1)
    	Y3 = .follow(n, Ltheta = .theta(Lhac$theta, ober.theta), I = First, type = type)
        v1 = .nextStage(n, d = dim, Y = Y3, Ltheta = Lhac$theta, type = type)
        colnames(v1) = Lhac$V1
        v2 = .simulate.hac(n, Lhac = Lhac$V2, First = Y3, ober.theta = Lhac$theta, type = type)
        v = cbind(v1, v2)}
    else if((class(Lhac$V1) != "character") & (class(Lhac$V2) != "character")){
    	Y4 = .follow(n, Ltheta = .theta(Lhac$theta, ober.theta), I = First, type = type)    
    	v1 = .simulate.hac(n, Lhac = Lhac$V1, First = Y4, ober.theta = Lhac$theta, type = type)
        v2 = .simulate.hac(n, Lhac = Lhac$V2, First = Y4, ober.theta = Lhac$theta, type = type)
        v = cbind(v1, v2)}
return(v)
}	

#-------------------------------------------------------------------------------------------------------------------------------

.rHAC = function(n, L){
hac = L$model
if((class(hac$V1) == "character") & (class(hac$V2) == "character")){
		Y1 = .initial(n, Ltheta = hac$theta, type = L$type)
		dim = length(hac$V1) + length(hac$V2)
		v  = .stayStage(n, d = dim, Y = Y1, Ltheta = hac$theta, type = L$type)
		colnames(v) = c(hac$V1, hac$V2)}
    else
if((class(hac$V1) != "character") & (class(hac$V2) == "character")){
		hac2 = hac$V1
		Y2 = .initial(n, Ltheta = hac$theta, type = L$type)
		dim = length(hac$V2)
  		v2 = .stayStage(n, d = dim, Y = Y2, Ltheta = hac$theta, type = L$type)
  		colnames(v2) = hac$V2
        v1 = .simulate.hac(n, Lhac = hac2, First = Y2, ober.theta = hac$theta, type = L$type)
		v  = cbind(v1, v2)}
    else
if((class(hac$V1) == "character") & (class(hac$V2) != "character")){
		hac3 = hac$V2
		dim = length(hac$V1)
		Y3 = .initial(n, Ltheta = hac$theta, type = L$type)
        v1 = .stayStage(n, d = dim, Y = Y3, Ltheta = hac$theta, type = L$type)
  		colnames(v1) = hac$V1        
        v2 = .simulate.hac(n, Lhac = hac3, First = Y3, ober.theta = hac$theta, type = L$type)
		v  = cbind(v1, v2)}
	else
if((class(hac$V1) != "character") & (class(hac$V2) != "character")){
		hac4 = hac$V1
		hac5 = hac$V2
		Y4 = .initial(n, hac$theta, L$type)
        v1 = .simulate.hac(n, Lhac = hac4, First = Y4, ober.theta = hac$theta, type = L$type)
        v2 = .simulate.hac(n, Lhac = hac5, First = Y4, ober.theta = hac$theta, type = L$type)
		v  = cbind(v1, v2)
        colnames(v) = c(colnames(v1), colnames(v2))}
	return(v)
}