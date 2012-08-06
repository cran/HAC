# copula_simulate.r ######################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .rAC						Samples from AC. (Internal function)
#  .f_gumbel				Samples from the inverse Laplace-Stietjes transfrom of a Gumbel copula. (Internal function)
#  .f_clayton				Samples from the inverse Laplace-Stietjes transfrom of a Clayton copula. (Internal function)
#  rHAC						Samples from HAC.
#  .theta      				Computes the ratio of two dependency parameters. (Internal function)
#  .initial  				Samples from the inverse Laplace-Stietjes transfrom for the initial node of HAC. (Internal function) 
#  .stayStage  			    Samples from the initial node of HAC. (Internal function)     
#  .fReject 				Samples from inverse Laplace-Stietjes transfrom for subsequent nodes, if type = HAC_CLAYTON. (Internal function)     
#  .follow					Samples the inverse Laplace-Stietjes transfrom of all successive nodes of nested AC. (Internal function)     
#  .simualte		  		The recursive sampling procedure. (Internal function)     
#  .rHAC 					Initializes the recursion. (Internal function)     
##########################################################################################################################

.rAC = function(n, theta = 1.5, dim = 2, type = AC_GUMBEL){
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
    tree = hac$tree
    type = hac$type
	if((type == AC_CLAYTON) | (type == AC_GUMBEL)){
		m = length(unlist(tree))
		res = .rAC(n, theta = tree[[m]], dim = (m-1), type = type)
        colnames(res) = unlist(tree)[-m]
	}else{
	if(hac$type == GAUSS){
        res = rCopula(n, normalCopula(tree[lower.tri(tree)], dim = NCOL(tree), dispstr = "un"))
	}else{
		res = .rHAC(n, tree, type)
	}}
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
		while(any((!Com) | (is.na(Com)))){
			fa = which((!Com) | (is.na(Com)))
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

.simulate = function(n, tree, First, ober.theta, type){
    dd = length(tree)
    m = matrix(, nrow = n)
    Y = .follow(n, Ltheta = .theta(tree[[dd]], ober.theta), I = First, type)
  
    select = sapply(tree[-dd], FUN = is.character)
    if(any(select)){
        this.node = which(select)
        v = as.matrix(.stayStage(n, d = length(this.node), Y = Y, Ltheta = tree[[dd]], type))
        colnames(v) = tree[this.node]
        m = cbind(m, v)
    }

    if(any(!select)){
    for(i in which(!select)){
       	v = .simulate(n, tree = tree[[i]], First = Y, ober.theta = tree[[dd]], type = type)
       	m = cbind(m, v)
	}}
	return(m)
}

#-------------------------------------------------------------------------------------------------------------------------------

.rHAC = function(n, tree, type){
    dd = length(tree)
    m = matrix(, nrow = n)
    Y = .initial(n, Ltheta = tree[[dd]], type = type)
    
    select = sapply(tree[-dd], FUN = is.character)
    if(any(select)){
        this.node = which(select)
        v = as.matrix(.stayStage(n, d = length(this.node), Y = Y, Ltheta = tree[[dd]], type = type))
        colnames(v) = tree[this.node]
        m = cbind(m, v)
    }
    
    if(any(!select)){
    for(i in which(!select)){
            v = .simulate(n, tree = tree[[i]], First = Y, ober.theta = tree[[dd]], type = type)
          	m = cbind(m, v)
	}}
	return(m[,-which(is.na(m[1,]))])
}
