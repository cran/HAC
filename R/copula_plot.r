# copula_plot.r ###########################################################################################################
# FUNCTION:               	DESCRIPTION:
#  plot.hac						Produces plots of the structure of HAC.
#  .subscribt					"Computes" the subsripts of the dependency parameter of the HAC. (Internal function)
#  .from.down					Calculates the coordinates of the nodes and variables of HAC. (Internal function)
#  .shift				 			Adjusts the coordinates of the variables. (Internal function)  
##########################################################################################################################


plot.hac = function(x, xlim = c(0, (dim(t(.get.leaves(tree)))[2]+1)), ylim = c(-0.5,(.coord(tree)$coord$y+0.5)), xlab = "", ylab = "", axes = FALSE, col = "black", fg = "black", bg = "white", lwd = 2, index = FALSE, theta = TRUE, h = 0.45, l = 1.2, circles = 0.25, digits = 2, ...){

.get.leaves = function(Ltree){
	if(class(Ltree) == "integer"){Ltree}
	else if( (class(Ltree$V1) != "character") & (class(Ltree$V2) != "character"))
		{c(.get.leaves(Ltree$V1), .get.leaves(Ltree$V2))}
	else if( (class(Ltree$V1) == "character") & (class(Ltree$V2) == "character"))
		{c(Ltree$V1, Ltree$V2)}
	else if( (class(Ltree$V1) == "character") & (class(Ltree$V2) != "character"))
		{c(Ltree$V1, .get.leaves(Ltree$V2))}
	else if( (class(Ltree$V1) != "character") & (class(Ltree$V2) == "character"))
		{c(.get.leaves(Ltree$V1), Ltree$V2)}
}

if((x$type == HAC_GUMBEL) | (x$type == HAC_ROTATED_GUMBEL) | (x$type == HAC_CLAYTON)){

tree = x$model

	A = ((FALSE == index) | (dim(t(.get.leaves(tree)))[2] > 7))
	B = (TRUE == index)
	C = c(A, B)
	s = 0.28 * (dim(t(.get.leaves(tree)))[2])
	
.coord = function(Ltree){
	if( (!("coord" %in% names(Ltree))) & (class(Ltree) != "character")){
        Ltree$V1 = .coord(Ltree$V1)
        Ltree$V2 = .coord(Ltree$V2)
		Ltree$coord <- .from.down(Ltree$V1$coord, Ltree$V2$coord)
	}else if(class(Ltree) == "character"){
		n = length(Ltree)
        Ltree = list(Ltree, coord = list(x = integer(n), y = integer(n)))
        for(i in 1 : n){
        	Ltree$coord$x[i] = which(Ltree[[1]][i] == .get.leaves(tree))
        }
    }
	Ltree
}

.circle = function(a, b, L){
	symbols(a, b, circles = circles, add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
	text(a, b, L[[1]])
	}
	
.lines = function(L){
	lines(c(L$coord$x,L$V1$coord$x), c(L$coord$y,L$V1$coord$y), lwd = lwd, col = col)
	lines(c(L$coord$x,L$V2$coord$x), c(L$coord$y,L$V2$coord$y), lwd = lwd, col = col)
}	

.rectangle = function(a, b, L, z, C, theta, type){
	if(theta == TRUE){
	switch(C, 
	A = {
		symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
		text(a, b, label = bquote(paste(theta == .(round(L$theta, digits = digits)))))},
	B = {
		symbols(a, b, rectangles = cbind(z, 0.45), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
		text(a, b, label = bquote(paste(theta[.(.subscribt(L))] == .(round(L$theta, digits = digits)))))})}
	else
	if(theta == FALSE){
	switch(C, 
	A = {
		symbols(a, b, rectangles = cbind(l, h), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
		text(a, b, label = bquote(paste(tau == .(round(theta2tau(L$theta, type), digits = digits)))))},
	B = {
		symbols(a, b, rectangles = cbind(z, 0.45), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
		text(a, b, label = bquote(paste(tau[.(.subscribt(L))] == .(round(theta2tau(L$theta, type), digits = digits)))))})
}}

#.plotlines = function(Ltree){
#	if(class(Ltree[[1]]) == "character"){.lines(Ltree)}
#		else{ if(class(Ltree$theta) == "numeric")
#			{.lines(Ltree)}	
#			.plotlines(Ltree$V1)
#			.plotlines(Ltree$V2)
#	}
#}

.plotlines = function(L){
	if((class(L$V1[[1]]) == "list") & (class(L$V2[[1]]) == "list")){
		.lines(L)
		.plotlines(L$V1)
		.plotlines(L$V2)
	}
	else
	if((class(L$V1[[1]]) == "list") & (class(L$V2[[1]]) == "character")){
		for(i in 1 : length(L$V2[[1]])){
			lines(c(L$coord$x, L$V2$coord$x[i]), c(L$coord$y, L$V2$coord$y[i]), lwd = lwd, col = col)
		}
		lines(c(L$coord$x,L$V1$coord$x), c(L$coord$y,L$V1$coord$y), lwd = lwd, col = col)
		.plotlines(L$V1)
	}
	else
	if((class(L$V1[[1]]) == "character") & (class(L$V2[[1]]) == "list")){
		for(i in 1 : length(L$V1[[1]])){
			lines(c(L$coord$x,L$V1$coord$x[i]), c(L$coord$y,L$V1$coord$y[i]), lwd = lwd, col = col)
		}
		lines(c(L$coord$x,L$V2$coord$x), c(L$coord$y,L$V2$coord$y), lwd = lwd, col = col)
		.plotlines(L$V2)
	}
	else
	if((class(L$V1[[1]]) == "character") & (class(L$V2[[1]]) == "character")){
		for(i in 1 : length(L$V1[[1]])){
			lines(c(L$coord$x,L$V1$coord$x[i]), c(L$coord$y,L$V1$coord$y[i]), lwd = lwd, col = col)
		}
		for(i in 1 : length(L$V2[[1]])){
			lines(c(L$coord$x,L$V2$coord$x[i]), c(L$coord$y,L$V2$coord$y[i]), lwd = lwd, col = col)
		}		
	}
}

.plottree = function(Ltree, L, l, C, theta, type){switch(C, 
		A = {
	if(class(Ltree[[1]]) == "character"){for(i in 1 : length(Ltree[[1]])){.circle(Ltree$coord$x[i], Ltree$coord$y[i], Ltree[[1]][i])}}
	else{if(class(Ltree$theta) == "numeric")
		{.rectangle(Ltree$coord$x, Ltree$coord$y, L, l, "A", theta, type)}
		.plottree(Ltree$V1, L$V1, l, "A", theta, type)
		.plottree(Ltree$V2, L$V2, l, "A", theta, type)
	}},  
		B = {
	if(class(Ltree[[1]]) == "character"){for(i in 1 : length(Ltree[[1]])){.circle(Ltree$coord$x[i], Ltree$coord$y[i], Ltree[[1]][i])}}
	else{if(class(Ltree$theta) == "numeric")
		{.rectangle(Ltree$coord$x, Ltree$coord$y, L, s, "B", theta, type)}
		.plottree(Ltree$V1, L$V1, s, "B", theta, type)
		.plottree(Ltree$V2, L$V2, s, "B", theta, type)}}
)}

if(A == TRUE){
	graphics::plot(x = 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, col = "white", ...)
	symbols(.coord(tree)$coord$x, .coord(tree)$coord$y, rectangles = cbind(l, h), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
	.plotlines(.shift(.coord(tree)))
	.plottree(.shift(.coord(tree)), tree, l, "A", theta, x$type)}
else{
	graphics::plot(x = 0, type = "p", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, col = "white", ...)
	symbols(.coord(tree)$coord$x, .coord(tree)$coord$y, rectangles = cbind(s ,0.25), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
	.plotlines(.shift(.coord(tree)))
	.plottree(.shift(.coord(tree)), tree, s, "B", theta, x$type)}
}
else
if((x$type == AC_GUMBEL) | (x$type == AC_CLAYTON)){
	model = x$model
	dim = x$model$dim
	the = x$model$theta
	tree = integer(dim)
	coord.circ = 1 : dim
	
	graphics::plot(1, col = "white", xlim = xlim + c(0.7, -0.7), ylim = ylim, xlab = xlab, ylab = ylab, axes = axes, ...)
	for(i in 1 : dim){
		lines(c(coord.circ[i], (max(coord.circ) + min(coord.circ))/2), c(0, 1), lwd = lwd, col = col)
		symbols(coord.circ[i], 0, circles = circles, add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
		text(coord.circ[i], 0, label = bquote(paste(V[.(coord.circ[i])])))
	}
	
	symbols((max(coord.circ) + min(coord.circ))/2, 1, rectangles = cbind(l, h), add = TRUE, inches = FALSE, bg = bg, lwd = lwd, fg = fg)
	if(theta == TRUE){
		text((max(coord.circ) + min(coord.circ))/2, 1, label = bquote(paste(theta == .(round(the, digits = digits)))))}
	else
	if(theta == FALSE){
		text((max(coord.circ) + min(coord.circ))/2, 1, label = bquote(paste(tau == .(round(theta2tau(the, x$type), digits = digits)))))}
}
}

#-------------------------------------------------------------------------------------------------------------------------------

.subscribt = function(Ltree){
    if((class(Ltree$V1) == "character") & (class(Ltree$V2) == "character"))
        {return(paste("(", paste(Ltree$V1, collapse = "."), ".", paste(Ltree$V2, collapse = "."), ")", sep = ""))}
    else 
    if((class(Ltree$V1) != "character") & (class(Ltree$V2) == "character"))
        {return(paste("(", .subscribt(Ltree$V1), ".(", paste(Ltree$V2, collapse = "."), "))", sep = ""))}
    else 
    if((class(Ltree$V1) == "character") & (class(Ltree$V2) != "character"))
        {return(paste("((", paste(Ltree$V1, collapse = "."), ").", .subscribt(Ltree$V2), ")", sep = ""))}
    else 
    if((class(Ltree$V1) != "character") & (class(Ltree$V2) != "character"))
        {return(paste("(", .subscribt(Ltree$V1), ".", .subscribt(Ltree$V2), ")", sep = ""))}
}

#-------------------------------------------------------------------------------------------------------------------------------

.from.down = function(L, R){
	result   = list()
	result$x = (min(L$x) + max(R$x))/2
	result$y =  max(L$y, R$y) + 1
	result
}

#-------------------------------------------------------------------------------------------------------------------------------

.shift = function(Ltree){
	n1 = length(Ltree$V1$coord$y)
	n2 = length(Ltree$V2$coord$y)
	Ltree$V1$coord$y = Ltree$V2$coord$y = (Ltree$coord$y-1)
	Ltree$V1$coord$y = rep(Ltree$V1$coord$y, n1)
	Ltree$V2$coord$y = rep(Ltree$V2$coord$y, n2)
	if(class(Ltree$V1[[1]]) != "character")Ltree$V1 = .shift(Ltree$V1)
	if(class(Ltree$V2[[1]]) != "character")Ltree$V2 = .shift(Ltree$V2)
	Ltree
	}