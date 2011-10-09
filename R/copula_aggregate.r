# copula_aggregate.r #####################################################################################################
# FUNCTION:               	DESCRIPTION:
#  .union					Tests, what kind of objects ((group of) variable(s) or / and node) are at the next nesting level. (Internal function)
#  .node.var           	Is called, if at the next nesting level are a node and a (group of) variable(s). (Internal function)
#  .var.node          		Is called, if at the next nesting level are a (group of) variable(s) and a node. (Internal function)
#  .node.node        		Is called, if at the next nesting level are two nodes. (Internal function)
#  as.hac				  	Transforms an object, which consists of 'list' object(s), into an hac object. (Internal function)
#  aggregate.hac         Aggregates an object of the class 'hac'.
# .method					Allows for the computation of mean, min and max if abs(theta_1 - theta_2) > epsilon             
##########################################################################################################################

.union = function(tree, epsilon, method){
 	if((class(tree$V1) == "list") & (class(tree$V2) == "character")){
 		tree = .node.var(tree = tree, epsilon = epsilon, method = method)}
 	else
 	if((class(tree$V1) == "character") & (class(tree$V2) == "list")){
 		tree = .var.node(tree = tree, epsilon = epsilon, method = method)}
 	else
 	if((class(tree$V1) == "character") & (class(tree$V2) == "character")){
 		tree$V1 = tree$V1
 		tree$V2 = tree$V2
 		tree}
 	else
 	if((class(tree$V1) == "list") & (class(tree$V2) == "list")){
 		tree = .node.node(tree = tree, epsilon = epsilon, method = method)}
 tree
}

#--------------------------------------------------------------------------------------------

.node.var = function(tree, epsilon, method){
 	if(abs(tree$theta - tree$V1$theta) < epsilon){
 		theta = .method(c(tree$theta, tree$V1$theta), method = method)
 			if((class(tree$V1$V1) == "character") & (class(tree$V1$V2) == "character")){
 				a = tree$V1$V1
 				b = tree$V1$V2
 				tree$V1 = c(a, b)
 				tree$theta = theta
 				tree
 			}
 		else
 			if((class(tree$V1$V1) == "list") & (class(tree$V1$V2) == "character")){
 				a = tree$V1$V2
 				b = tree$V2
 				tree$V2 = c(a, b)
 				tree$theta = theta
 				tree$V1 = tree$V1$V1
				tree = .union(tree, epsilon = epsilon, method = method)
 			}
 		else
 			if((class(tree$V1$V1) == "character") & (class(tree$V1$V2) == "list")){
 				a = tree$V2
 				b = tree$V1$V1
 				tree$V2 = c(a, b)
 				tree$theta = theta
 				tree$V1 = tree$V1$V2
 				tree = .union(tree, epsilon = epsilon, method = method)	
 			}
 		else
 			if((class(tree$V1$V1) == "list") & (class(tree$V1$V2) == "list")){
 				tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = tree$V2, theta = tree$theta)
 			}
 	}
 else
 	if(abs(tree$theta - tree$V1$theta) >= epsilon){
 		tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = tree$V2, theta = tree$theta)
 	}
}
 
#--------------------------------------------------------------------------------------------
 
.var.node = function(tree, epsilon, method){
 	if(abs(tree$theta - tree$V2$theta) < epsilon){
 		theta = .method(c(tree$theta, tree$V2$theta), method = method)
			if((class(tree$V2$V1) == "character") & (class(tree$V2$V2) == "character")){
 				a = tree$V2$V1
 				b = tree$V2$V2
 				tree$V2 = c(a, b)
 				tree$theta = theta
 				tree
 			}
 			else
 			if((class(tree$V2$V1) == "list") & (class(tree$V2$V2) == "character")){
 				a = tree$V2$V2
 				b = tree$V1
 				tree$V1 = c(a, b)
 				tree$theta = theta
 				tree$V2 = tree$V2$V1
 				tree = .union(tree, epsilon = epsilon, method = method)
 			}
 			else
 			if((class(tree$V2$V1) == "character") & (class(tree$V2$V2) == "list")){
 				a = tree$V1
 				b = tree$V2$V1
 				tree$V1 = c(a, b)
 				tree$theta = theta
 				tree$V2 = tree$V2$V2
 				tree = .union(tree, epsilon = epsilon, method = method)
 			}
 			else
 			if((class(tree$V2$V1) == "list") & (class(tree$V2$V2) == "list")){
				tree = list(V1 = tree$V1, V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
 			}
 	}
 else
 	if(abs(tree$theta - tree$V2$theta) >= epsilon){
 		tree = list(V1 = tree$V1, V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
 	}
}
 
#--------------------------------------------------------------------------------------------
 
.node.node = function(tree, epsilon, method){
 	if((abs(tree$theta - tree$V1$theta) < epsilon) & (abs(tree$theta - tree$V2$theta) < epsilon)){
  		theta = .method(c(.method(c(tree$theta, tree$V1$theta), method = method), .method(c(tree$theta, tree$V2$theta), method = method)), method = method)
  			if((class(tree$V1$V1) == "character") & (class(tree$V1$V2) == "character") & (class(tree$V2$V1) == "character") & (class(tree$V2$V2) == "character")){
  				tree$theta = theta
  				a = tree$V1$V1
  				b = tree$V1$V2
  				c = tree$V2$V1
  				d = tree$V2$V2
  				tree$V1 = c(a, b)
  				tree$V2 = c(c, d)
  				tree}
  		else
  			if(((class(tree$V1$V1) == "list") | (class(tree$V1$V2) == "list")) & ((class(tree$V2$V1) == "character") & (class(tree$V2$V2) == "character"))){
  				tree$theta = theta
  				c = tree$V2$V1
  				d = tree$V2$V2
  				tree$V2 = c(c, d)
  				tree = .union(tree, epsilon = epsilon, method = method)}
  		else
  			if(((class(tree$V1$V1) == "character") & (class(tree$V1$V2) == "character")) & ((class(tree$V2$V1) == "list") | (class(tree$V2$V2) == "list"))){
  				tree$theta = theta
  				a = tree$V1$V1
  				b = tree$V1$V2
  				tree$V1 = c(a, b)
  				tree = .union(tree, epsilon = epsilon, method = method)}
  		else
  			if(((class(tree$V1$V1) == "list") | (class(tree$V1$V2) == "list")) & ((class(tree$V2$V1) == "list") | (class(tree$V2$V2) == "list"))){
  				tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
  			}}
  	else
  		if((abs(tree$theta - tree$V1$theta) >= epsilon) & (abs(tree$theta - tree$V2$theta) < epsilon)){
 			theta = .method(c(tree$theta, tree$V2$theta), method = method)
  			if((class(tree$V2$V1) == "character") & (class(tree$V2$V2) == "character")){
  				tree$theta = theta
  				a = tree$V2$V1
				b = tree$V2$V2
  				tree$V2 = c(a, b)
  				tree = .union(tree, epsilon = epsilon, method = method)
  			}
  		else
 			if((class(tree$V2$V1) != "character") | (class(tree$V2$V2) != "character")){
  				tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
  			}
  			}
  	else
  		if((abs(tree$theta - tree$V1$theta) < epsilon) & (abs(tree$theta - tree$V2$theta) >= epsilon)){
 			theta = .method(c(tree$theta, tree$V1$theta), method = method)
  			if((class(tree$V1$V1) == "character") & (class(tree$V1$V2) == "character")){
  				a = tree$V1$V1
  				b = tree$V1$V2
 				tree$theta = theta
  				tree$V1 = c(a, b)
  				tree = .union(tree, epsilon = epsilon, method = method)
  			}
  		else
  			if((class(tree$V1$V1) != "character") | (class(tree$V1$V2) != "character")){
  				tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
 				}		
  			}
  	else
 		if((abs(tree$theta - tree$V1$theta) >= epsilon) & (abs(tree$theta - tree$V2$theta) >= epsilon)){
			tree = list(V1 = .union(tree$V1, epsilon = epsilon, method = method), V2 = .union(tree$V2, epsilon = epsilon, method = method), theta = tree$theta)
 		}
 }
 
#--------------------------------------------------------------------------------------------
 
as.hac = function(type = HAC_GUMBEL, tree){
 	object = list(type = type, model = tree)
 	class(object) = "hac"
 	object
}

#--------------------------------------------------------------------------------------------

aggregate.hac = function(x, epsilon = 0.01, method = "mean", ...){
 	tree = .union(x$model, epsilon = epsilon, method = method)
 	if((class(tree$V1) == "character") & (class(tree$V2) == "character")){
 		if(x$type == HAC_GUMBEL){
 			hac(type = AC_GUMBEL, X = tree$theta, dim = (length(tree$V1) + length(tree$V2)))}
		else
		if(x$type == HAC_ROTATED_GUMBEL){
			hac(type = HAC_ROTATED_GUMBEL, X = tree$theta, dim = (length(tree$V1) + length(tree$V2)))}
		else
		if(x$type == HAC_CLAYTON){
			hac(type = AC_CLAYTON, X = tree$theta, dim = (length(tree$V1) + length(tree$V2)))}				
	}
else{
	as.hac(type = x$type, tree = tree)
	}	
}

#--------------------------------------------------------------------------------------------

.method = function(x, method){
	if(method == "mean"){mean(x)}
	else
	if(method == "min"){min(x)}
	else
	if(method == "max"){max(x)}
}