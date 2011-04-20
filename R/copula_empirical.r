# copula_empirical.r #####################################################################################################
# FUNCTION:               	DESCRIPTION:
#  emp.copula				Computes the emprical copula for a given grid u and sample x.
#  emp.copula.self          Computes the emprical copula for a given sample x.       
##########################################################################################################################

emp.copula = function(u, x, proc = "M", sort = "none", na.rm = FALSE){
	
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
	x.help = matrix((rowSums(is.na(x < 0)) + rowSums(is.na(x > 1)) + rowSums(is.na(x))), nrow = n)
	if(sum(x.help) != 0){
	x = x[-which(x.help > 0), ]
	for(i in 1 : n){
	if(x.help[i,] > 0){warning(paste("Row", i,"of x was deleted, since at least one element of the row is >", 1,", <", 0,"or is.na(x[", i,", ]) == TRUE."))}}}}

	n = dim(x)[1]
	d = dim(x)[2]
	
	if(na.rm == TRUE){
	u.help = matrix((rowSums(is.na(u < 0)) + rowSums(is.na(u > 1)) + rowSums(is.na(u))), nrow = nn)
	if(sum(u.help) != 0){
	u = u[-which(u.help > 0), ]
	for(i in 1 : nn){
	if(u.help[i,] > 0){warning(paste("Row", i,"of u was deleted, since at least one element of the row is >", 1,", <", 0,"or is.na(u[", i,", ]) == TRUE."))}}}}
	
	nn = dim(u)[1]
	dd = dim(u)[2]
			
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

emp.copula.self = function(x, proc = "M", sort = "none", na.rm = FALSE){switch(proc, M = emp.copula(x, x, "M", sort = sort, na.rm = na.rm), A = emp.copula(x, x, "A", sort = sort, na.rm = na.rm))}
