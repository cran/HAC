# copula_empirical.r #####################################################################################################
# FUNCTION:               	DESCRIPTION:
#  emp.copula				Computes the emprical copula for a given grid u and sample x.
#  emp.copula.self          Computes the emprical copula for a given sample x.       
##########################################################################################################################

emp.copula = function(u, x, proc = "M", sort = "none", na.rm = FALSE, max.min = TRUE){
	
	if((class(u) != "matrix") & (class(u) == "numeric")){
		 u = t(u)}
	
	 n = dim(x)[1]
	 d = dim(x)[2]
	nn = dim(u)[1]
	dd = dim(u)[2]
	
	if((dd != d) & (nn == d)){
		nn = dd
		 u = t(u)}
	
	if(na.rm == TRUE){
		if(any(is.na(x))){
			x.help = matrix(rowSums(is.na(x)), nrow = n)
			rm = which(x.help > 0)
			x = x[-rm, ]
			n = dim(x)[1]
			warning(paste("The following rows of x were removed, since at least one element of them is NA:"), paste(rm, collapse = ","))}}
	
	if(na.rm == TRUE){
		if(any(is.na(u))){
			u.help = matrix(rowSums(is.na(u)), nrow = nn)
			rm = which(u.help > 0)
			u = u[-rm, ]
			nn = dim(u)[1]
			warning(paste("The following rows of u were removed, since at least one element of them is NA:"), paste(rm, collapse = ","))}}

	if(max.min == TRUE){
		x[which(x > 1)] = 1
		x[which(x < 0)] = 0
		u[which(u > 1)] = 1
		u[which(u < 0)] = 0}
	
	emp = function(u, x, proc){
		if(proc == "M"){
		Compare = matrix(t(matrix(rep(x, nn), ncol = nn * d)), ncol = d, byrow = TRUE)
		Values = matrix(rep(t(x), nn), n * nn, d, byrow = TRUE)
		
		1 / n * colSums(matrix((rowSums(Values <= Compare) == d), ncol = nn)) 
	}
	
	else if(proc == "A"){
		
		embedded = function(u){
			compare.help = function(X){(X <= u)}
				Compare = function(X){apply(X, 1, compare.help)}
			1 / n * cumsum(apply(t(Compare(x)), 1, prod))
		}
		
		apply(u, 1, embedded)[n,]
	}}

if(sort == "none"){
	emp(u = u, x = x, proc = proc)}
else
if(sort == "asc"){
	sort(emp(u = u, x = x, proc = proc))}
else
if(sort == "desc"){
	sort(emp(u = u, x = x, proc = proc), decreasing = TRUE)}
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

emp.copula.self = function(x, proc = "M", sort = "none", na.rm = FALSE, max.min = TRUE){switch(proc, M = emp.copula(x, x, "M", sort = sort, na.rm = na.rm, max.min = max.min), A = emp.copula(x, x, "A", sort = sort, na.rm = na.rm, max.min =
max.min))}
