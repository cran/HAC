# copula_empirical.r #####################################################################################################
# FUNCTION:               	DESCRIPTION:
#  emp.copula				Returns the values of the emprical copula for a given grid u and sample x.
#  emp.copula.self        	Returns the values of the emprical copula for a given sample x.  
#  .emp   						Computes the values of the empirical copula. (Internal function)
##########################################################################################################################

emp.copula = function(u, x, proc = "M", sort = "none", margins = NULL, na.rm = FALSE, max.min = TRUE){
    
	if((class(u) != "matrix") & (class(u) == "numeric")){
		 u = t(u)}
	
	 n = dim(x)[1]
	 d = dim(x)[2]
	nn = dim(u)[1]
	dd = dim(u)[2]
	
	if((dd != d) & (nn == d)){
		nn = dd
		 u = t(u)}
	
	x = .margins(x, margins)
	
	if(na.rm == TRUE){
		x = .na.rm(X = x)
		u = .na.rm(X = u)}

	if(max.min == TRUE){
		x = .max.min(X = x)
		u = .max.min(X = u)}
	
if(sort == "none"){
	.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d)}
else
if(sort == "asc"){
	sort(.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d))}
else
if(sort == "desc"){
	sort(.emp(u = u, x = x, proc = proc, n = n, nn = nn, d = d), decreasing = TRUE)}
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

emp.copula.self = function(x, proc = "M", sort = "none", margins = NULL, na.rm = FALSE, max.min = TRUE){switch(proc, M = emp.copula(x, x, "M", sort = sort, margins = margins, na.rm = na.rm, max.min = max.min), A = emp.copula(x, x, "A", sort = sort, margins = margins, na.rm = na.rm, max.min =
max.min))}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

.emp = function(u, x, proc, n, nn, d){
	if(proc == "M"){
		Compare = matrix(t(matrix(rep(x, nn), ncol = nn * d)), ncol = d, byrow = TRUE)
		Values = matrix(rep(t(x), nn), n * nn, d, byrow = TRUE)
		
		1 / n * colSums(matrix((rowSums(Values <= Compare) == d), ncol = nn)) 
	}else{ if(proc == "A"){
	
	embedded = function(u){
		compare.help = function(X){(X <= u)}
		Compare = function(X){apply(X, 1, compare.help)}
		1 / n * cumsum(apply(t(Compare(x)), 1, prod))
	}
		apply(u, 1, embedded)[n,]
	}}
}