# copula_constructor.r ###################################################################################################
# FUNCTION:               	DESCRIPTION:
#  columns.hac				Calculates the maximal number of columns for the matrix X (argument of the function hac).
#  hac        				Constructs arbitrary, (binary) Hierarchical Archimedean Copulae.
#  .app1       				Returns the values of the matrix X. (Internal function)
#  .app2					Converts the matrix X in a (copula)-tree. (Internal function)
#  .max.model   			Constructs the biggest possible binary copulae frame for a given number of rows of X. (Internal function)
#  .is.wholenumber			Tests, whether its argument is an integer value. (Internal function)
#  check.var				Tests, whether each path of the copula ends with a variable.
#  check.par				Tests, whether the parameter are correctly specified.
#  hac.full					Constructs fully nested Archimedean Copulae.
#  print.hac     			Determines how an object of the class 'hac' is printed.
##########################################################################################################################

columns.hac = function(n){
	if(n < 2){return(warning("The argument has to be bigger than 1"))}
	if(n == 2){
		2}
	else{
	if(n > 2){
	2 * columns.hac(n - 1)}}
}

#------------------------------------------------------------------------------------------------------------------------

hac = function(type = HAC_GUMBEL, X, dim = NULL){
	
	if((NROW(X)==1) & (NCOL(X)==1) & (is.null(dim)==FALSE)){
		if(type == HAC_GUMBEL){type = AC_GUMBEL}
			else
			if(type == HAC_CLAYTON){type = AC_CLAYTON}
			
		CopulaNew = list(dim = 0, theta = 0)
		CopulaNew$dim = dim
		CopulaNew$theta = X
	
		result = list()
		result = list(type = type, model = CopulaNew)
		class(result) = "hac"
		result
	}
	else
	if((type == HAC_GUMBEL) | (type == HAC_ROTATED_GUMBEL) | (type == HAC_CLAYTON)){
	
	n = dim(X)[1]
	h = dim(X)[2]

matrix = matrix(c(X, matrix(0, nrow = n, ncol = (columns.hac(n) - h))), ncol = columns.hac(n))

	d = dim(matrix)[2]

for(i in 2 : d){
		if(matrix[1, i][[1]] != 0)
		return(warning(paste("The matrix is not correctly specified. Element [ 1,", i, "] of the matrix needs to be equal to zero or NULL. Check the appropriate specification of the matrix and apply the function hac again to obtain the desired object.")))
		}
if(n > 2){
for(j in 2 : (n - 1)){
	for(i in  (2 * (j - 1) + 1) : d){
		if(matrix[j, i][[1]] != 0)
		return(warning(paste("The matrix is not correctly specified. Element [", j, ",", i, "] of the matrix needs to be equal to zero or NULL. Check the appropriate specification of the matrix and apply the function hac again to obtain the desired object.")))
	}}}
if(n > 1){
for(j in 1 : (dim(X)[1] - 1)){
	for(i in  1 : (dim(X)[2] / 2)){
		if((X[j, i][[1]] > 0) & (class(X[j, i][[1]]) == "numeric") & (X[(j + 1), (2 * i - 1)][[1]] == 0))
		return(warning(paste("The matrix is not correctly specified. Element [", (j + 1), ",", (2 * i - 1), "] of the matrix needs to be !=", 0, "Check the appropriate specification of the matrix and apply the function hac again to obtain the desired object.")))
	}
}}
if(n > 1){
for(j in 1 : (dim(X)[1] - 1)){
	for(i in  1 : (dim(X)[2] / 2)){
		if((X[j, i][[1]] > 0) & (class(X[j, i][[1]]) == "numeric") & (X[(j + 1), (2 * i)][[1]] == 0))
		return(warning(paste("The matrix is not correctly specified. Element [", (j + 1), ",", (2 * i), "] of the matrix needs to be !=", 0, "Check the appropriate specification of the matrix and apply the function hac again to obtain the desired object.")))
	}}
}

	cmodel = list()
	cmodel = list(type = type, model = .app2(.max.model(n), 1, 1, matrix))
	check.var(cmodel)
	check.par(cmodel)
	class(cmodel) = "hac"
	cmodel}
}

#------------------------------------------------------------------------------------------------------------------------

.app1 = function(j, i, X){
if((.is.wholenumber(i / 2) == FALSE) & (class(X[j, i][[1]]) == "formula") & (X[j, i][[1]] != 0)){
	attr(terms(X[j, i][[1]]), "term.labels")}
else
if((.is.wholenumber(i / 2) == TRUE) & (class(X[j, i][[1]]) == "formula") & (X[j, i][[1]] != 0)){
	attr(terms(X[j, i][[1]]), "term.labels")}
else
if((.is.wholenumber(i / 2) == FALSE) & (class(X[j, i][[1]]) != "formula") & (X[j, i][[1]] != 0)){
	X[j, i][[1]]}
else
if((.is.wholenumber(i / 2) == TRUE) & (class(X[j, i][[1]]) != "formula") & (X[j, i][[1]] != 0)){
	X[j, i][[1]]}
}

#------------------------------------------------------------------------------------------------------------------------

.max.model = function(z){
	.tree = function(Copula){
		CopulaNew = list(V1 = Copula, V2 = Copula, theta = 0)
		CopulaNew	
		}

	Cop = list(V1 = list(), V2 = list(), theta = 0)

	if(z <= 2)Cop
	else if(z == 3).tree(Cop)
	else if(z >  3).tree(.max.model(z - 1))
}

#------------------------------------------------------------------------------------------------------------------------

.is.wholenumber = function(X, tol = .Machine$double.eps^0.5){abs(X - round(X)) < tol}

#------------------------------------------------------------------------------------------------------------------------

.app2 = function(L, j, i, X){
	mapping = function(L, j, i, X){
		if((class(.app1(j, i, X)) == "character")){
			L = .app1(j, i, X)}
		else{
		if(class(.app1(j, i, X)) == "numeric"){
				L$theta = .app1(j, i, X)
			if(.is.wholenumber(i / 2) == FALSE){
				L$V1 = mapping(L$V1, j + 1, (2 * i - 1), X)
				L$V2 = mapping(L$V2, j + 1, (2 * i - 0), X)}
			if(.is.wholenumber(i / 2) == TRUE){
				L$V1 = mapping(L$V1, j + 1, (2 * i - 1), X)
				L$V2 = mapping(L$V2, j + 1, (2 * i - 0), X)}
			}}
		L
	}
mapping(L, j, i, X)
}

#------------------------------------------------------------------------------------------------------------------------

check.var = function(x){
	.check.var = function(L){
		if ((class(L) == "list") & (length(L) == 0)){
			return(warning("The model is not correctly specified. A path ends with an empty object."))}
		else{
		if((class(L) == "list") & (length(L) != 0)){
			.check.var(L$V1)
			.check.var(L$V2)}}
	}
	L = x$model
	.check.var(L = L)
}

#------------------------------------------------------------------------------------------------------------------------

check.par = function(x){
	if((x$type == AC_GUMBEL) || (x$type == HAC_GUMBEL) || (x$type == HAC_ROTATED_GUMBEL)){
		ober.theta = 1}
	else
	if((x$type = AC_CLAYTON) || (x$type = HAC_CLAYTON)){
		ober.theta = 0}
	
check = function(L, theta){
		if(L$theta < theta){
			return(warning("The dependency parameter of a deeper nested AC should have a higher value."))}
		else{
		if((class(L$V1) == "list") && (class(L$V2) == "list")){
			check(L$V1, theta = L$theta)
			check(L$V2, theta = L$theta)}
		else
		if((class(L$V1) == "list") && (class(L$V2) != "list")){
				check(L$V1, theta = L$theta)}
		else
		if((class(L$V1) != "list") && (class(L$V2) == "list")){
				check(L$V2, theta = L$theta)}
		}
	}
	L = x$model
	check(L = L, theta = ober.theta)
}

#------------------------------------------------------------------------------------------------------------------------

hac.full = function(type = HAC_GUMBEL, y, theta){
	if(length(y) != (length(theta) + 1)){
	return(warning("The input arguments does not fit to a fully nested HAC"))}
	
	M = matrix(c(matrix(c(theta, y[max(length(y))]), ncol = 1), 0, y[-max(length(y))]), ncol = 2)
	hac(type = type, X = M)
}

#------------------------------------------------------------------------------------------------------------------------
	
print.hac = function(x, digits = 2, ...){
		cat("Class: hac", "\n", ...)
	if((x$type == 0) | (x$type == 1) | (x$type == 2)){
 		cat("Generator: Gumbel", "\n", ...)
 		cat(tree2str(x, digits = digits),  "\n", ...)
 		}
	if((x$type == 3) | (x$type == 4)){
 		cat("Generator: Clayton", "\n", ...)
 		cat(tree2str(x, digits = digits),  "\n", ...)
 		}
 }