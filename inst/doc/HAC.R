### R code from vignette source 'HAC.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("HAC")


###################################################
### code chunk number 2: fourfully
###################################################
Obj1 = hac.full(type = HAC_GUMBEL, y = c("u4", "u3", "u2", "u1"), theta = c(2, 3, 4))
par(mai = c(0, 0, 0, 0))
plot(Obj1, index = TRUE, l = 1.6)


###################################################
### code chunk number 3: fourpartially
###################################################
Obj2 = hac(tree = list(list("u4", "u3", 3), list("u1", "u2", 4), 2))
par(mai = c(0, 0, 0, 0))
plot(Obj2, index = TRUE, l = 1.6)


###################################################
### code chunk number 4: HAC.Rnw:280-284
###################################################
library("HAC")
data("finData")
system.time(result <- estimate.copula(finData, margins = "edf"))
result


###################################################
### code chunk number 5: Scatter1
###################################################
par(mai = c(0, 0, 0, 0))
pairs(finData, pch = 20)


###################################################
### code chunk number 6: HAC.Rnw:360-361
###################################################
names(formals(estimate.copula))


###################################################
### code chunk number 7: result-agg
###################################################
    result.agg = estimate.copula(finData, margins = "edf", epsilon = 0.3)
    par(mai = c(0, 0, 0, 0))
    plot(result.agg, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 8: result
###################################################
    par(mai = c(0, 0, 0, 0))
    plot(result, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 9: HAC.Rnw:392-395 (eval = FALSE)
###################################################
## result.agg = estimate.copula(sample, margins = "edf", epsilon = 0.3)
## plot(result, circles = 0.3, index = TRUE, l = 1.7)
## plot(result.agg, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 10: HAC.Rnw:436-440
###################################################
G.cop = hac.full(type  = HAC_GUMBEL,
                 y     = c("X4", "X3", "X2", "X1"),
                 theta = c(1.1, 1.8, 2.5))
G.cop


###################################################
### code chunk number 11: HAC.Rnw:447-448
###################################################
hac(tree = list("X1", "X2", "X3", "X4", 2))


###################################################
### code chunk number 12: HAC.Rnw:451-452
###################################################
hac(tree = list(list("X1", "X2", 2.5), "X3", "X4", 1.5))


###################################################
### code chunk number 13: HAC.Rnw:456-460
###################################################
HAC = hac(tree = list(list("Y1", list("Z3", "Z4", 3), "Y2", 2.5),
                      list("Z1", "Z2", 2), list("X1", "X2", 2.4),
                      "X3", "X4", 1.5))
HAC


###################################################
### code chunk number 14: HAC
###################################################
par(mai = c(0, 0, 0, 0))
plot(HAC, cex = 0.8, circles = 0.35)


###################################################
### code chunk number 15: HAC.Rnw:473-474
###################################################
plot(HAC, cex = 0.8, circles = 0.35)


###################################################
### code chunk number 16: HAC.Rnw:487-488
###################################################
names(formals(plot.hac))


###################################################
### code chunk number 17: Scatter2
###################################################
set.seed(1)
sim.data = rHAC(500, G.cop)
par(mai = c(0, 0, 0, 0))
pairs(sim.data, pch = 20)


###################################################
### code chunk number 18: HAC.Rnw:526-528 (eval = FALSE)
###################################################
## sim.data = rHAC(500, G.cop)
## pairs(sim.data, pch = 20)


###################################################
### code chunk number 19: HAC.Rnw:544-545
###################################################
probs = pHAC(X = sim.data, hac = G.cop)


###################################################
### code chunk number 20: pp
###################################################
probs.emp = emp.copula.self(sim.data, proc = "M")
plot(probs, probs.emp, pch = 20, xlab = "True Probabilities", ylab = "Empirical Probabilites", asp = 1)
grid(lwd = 2)
points(probs, probs.emp, pch = 20)
lines(c(0,1), c(0,1), col = "red3", lwd = 2)


###################################################
### code chunk number 21: HAC.Rnw:572-573 (eval = FALSE)
###################################################
## probs.emp = emp.copula.self(sim.data, proc = "M")


###################################################
### code chunk number 22: HAC.Rnw:578-582 (eval = FALSE)
###################################################
## emp.copula(u, x, proc = "M", sort = "none", margins = NULL,
##            na.rm = FALSE, ...)
## emp.copula.self(x, proc = "M", sort = "none", margins = NULL,
##                 na.rm = FALSE, ...)


###################################################
### code chunk number 23: speed2 (eval = FALSE)
###################################################
## 
## ####
## #Figure6---------------------------------------------------------------------------------------------------------------------- #Note, that the result depends on the CPU you use
## 
## d = 5
## s = seq(550, 6275, by=25)
## set.seed(1)
## data = matrix(runif(max(s)*d), ncol = d, nrow = max(s))
## t = matrix(, nrow = length(s)-1, ncol = 2)
## for(i in 2:length(s)){
##     t1 = Sys.time()
##     a1 = emp.copula.self(as.matrix(data[1:s[i],]), proc = "M")
##     t2 = Sys.time();
##     a2 = emp.copula.self(as.matrix(data[1:s[i],]), proc = "A")
##     t3 = Sys.time()
##     t[i-1, 1] = difftime(t2, t1)
##     t[i-1, 2] = difftime(t3, t2)
## }
## 
## times_plot2 = t
##     ## do not remove # in front of pdf() and dev.off()
##     pdf("HAC-speed2.pdf", width = 8, height = 6)
## plot(s[2:190], t[1:189,1], ylim = c(t[1,1], t[206,2]), xlim = c(s[2], s[206]), xlab = "observations", ylab = "seconds",
## type = "l", lwd = 2, log = "xy", asp = 1)
## points(s[2:230], t[,2], type = "l", lty = 2, lwd = 2)
## grid()
##     dev.off()


###################################################
### code chunk number 24: HAC.Rnw:637-984 (eval = FALSE)
###################################################
## 
## #####Table 4---------------------------------------------------------------------------------------------
## 
## library("HAC")
## 
## ####
## # Construct the 3-dimensional hac objects for Gumbel and Clayton
##     hac3_G = hac.full(type = HAC_GUMBEL, y = c("X1", "X2", "X3"), theta = c(tau2theta(c(1/3, 2/3))))
##     hac3_C = hac.full(type = HAC_CLAYTON, y = c("X1", "X2", "X3"), theta = c(tau2theta(c(1/3, 2/3), type = HAC_CLAYTON)))
## 
## ####
## # Define the number of estimates times = 1000, the sample size n = 250 and prepare matrix to save the results
##     times = 1000; n = 250
##     results3_G = matrix(, nrow = times, ncol = 3)
##     results3_C = matrix(, nrow = times, ncol = 3)
##     colnames(results3_G) = c("theta1", "theta2", "structure")
##     colnames(results3_C) = c("theta1", "theta2", "structure")
## 
## ####
## # Simulation and HAC-estimation is repeated 1000 times within the loop
## for(i in 1:times){
##     data_hac3_G = rHAC(n, hac3_G)
##     data_hac3_C = rHAC(n, hac3_C)
## 
##     a = estimate.copula(data_hac3_G, type = HAC_GUMBEL, method = ML)
##     results3_G[i, 1:2] = get.params(a, sort = TRUE, decreasing = TRUE)
## 
##     # Check whether the structure is correct
##     results3_G[i, 3] =
##     if(class(a$tree[[1]])=="character"){if(a$tree[[1]]=="X1"){TRUE}else{FALSE}}else{if(a$tree[[2]]=="X1"){TRUE}else{FALSE}}
## 
##     a = estimate.copula(data_hac3_C, type = HAC_CLAYTON, method = ML)
##     results3_C[i, 1:2] = get.params(a, sort = TRUE, decreasing = TRUE)
## 
##     # Check whether the structure is correct
##     results3_C[i, 3] =
##     if(class(a$tree[[1]])=="character"){if(a$tree[[1]]=="X1"){TRUE}else{FALSE}}else{if(a$tree[[2]]=="X1"){TRUE}else{FALSE}}
## }
## 
## ###
## # Print the results for the 3-dimensional models
## summary(results3_G)
## apply(results3_G, 2, sd)
## 
## summary(results3_C)
## apply(results3_C, 2, sd)
## 
## ####
## # Construct the 5-dimensional hac objects for Gumbel and Clayton
## hac5_G = hac.full(type = HAC_GUMBEL, y = c("X1", "X2", "X3", "X4", "X5"), theta = c(tau2theta(c(1/9, 3/9, 5/9, 7/9))))
## hac5_C = hac.full(type = HAC_CLAYTON, y = c("X1", "X2", "X3", "X4", "X5"), theta = c(tau2theta(c(1/9, 3/9, 5/9, 7/9), type = HAC_CLAYTON)))
## 
## ####
## # Define all structures, which correspond to the same copula as hac5_G and hac5_C respectively and save them as string, i.e., struc1, struc2,...
## 
## hac5_G1 = hac(tree = list("X1", list("X2", list(list("X4", "X5", 4), "X3", 3) ,2.2), 1))
## struc1 = tree2str(hac5_G1, theta=FALSE)
## hac5_G2 = hac(tree = list("X1", list("X2", list(list("X5", "X4", 4), "X3", 3) ,2.2), 1))
## struc2 = tree2str(hac5_G2, theta=FALSE)
## hac5_G3 = hac(tree = list("X1", list("X2", list("X3", list("X4", "X5", 4), 3) ,2.2), 1))
## struc3 = tree2str(hac5_G3, theta=FALSE)
## hac5_G4 = hac(tree = list("X1", list("X2", list("X3", list("X5", "X4", 4), 3) ,2.2), 1))
## struc4 = tree2str(hac5_G4, theta=FALSE)
## hac5_G5 = hac(tree = list("X1", list(list(list("X4", "X5", 4), "X3", 3), "X2",2.2), 1))
## struc5 = tree2str(hac5_G5, theta=FALSE)
## hac5_G6 = hac(tree = list("X1", list(list(list("X5", "X4", 4), "X3", 3), "X2",2.2), 1))
## struc6 = tree2str(hac5_G6, theta=FALSE)
## hac5_G7 = hac(tree = list("X1", list(list("X3", list("X4", "X5", 4), 3), "X2",2.2), 1))
## struc7 = tree2str(hac5_G7, theta=FALSE)
## hac5_G8 = hac(tree = list("X1", list(list("X3", list("X5", "X4", 4), 3), "X2",2.2), 1))
## struc8 = tree2str(hac5_G8, theta=FALSE)
## hac5_G9 = hac(tree = list(list("X2", list(list("X4", "X5", 4), "X3", 3) ,2.2), "X1", 1))
## struc9 = tree2str(hac5_G9, theta=FALSE)
## hac5_G10 = hac(tree = list(list("X2", list(list("X5", "X4", 4), "X3", 3) ,2.2), "X1", 1))
## struc10 = tree2str(hac5_G10, theta=FALSE)
## hac5_G11 = hac(tree = list(list("X2", list("X3", list("X4", "X5", 4), 3) ,2.2), "X1", 1))
## struc11 = tree2str(hac5_G11, theta=FALSE)
## hac5_G12 = hac(tree = list(list("X2", list("X3",list("X5", "X4", 4), 3) ,2.2), "X1", 1))
## struc12 = tree2str(hac5_G12, theta=FALSE)
## hac5_G13 = hac(tree = list(list(list(list("X4", "X5", 4), "X3", 3), "X2",2.2), "X1", 1))
## struc13 = tree2str(hac5_G13, theta=FALSE)
## hac5_G14 = hac(tree = list(list(list(list("X5", "X4", 4), "X3", 3), "X2",2.2), "X1", 1))
## struc14 = tree2str(hac5_G14, theta=FALSE)
## hac5_G15 = hac(tree = list(list(list("X3", list("X4", "X5", 4), 3), "X2",2.2), "X1", 1))
## struc15 = tree2str(hac5_G15, theta=FALSE)
## hac5_G16 = hac(tree = list(list(list("X3", list("X5", "X4", 4), 3), "X2",2.2), "X1", 1))
## struc16 = tree2str(hac5_G16, theta=FALSE)
## 
## ####
## # Define the initial matrices for the results
## results5_G = matrix(0, nrow = 1, ncol = 5)
## results5_C = matrix(0, nrow = 1, ncol = 5)
## colnames(results5_G) = c("theta1", "theta2", "theta3", "theta4", "structure")
## colnames(results5_C) = c("theta1", "theta2", "theta3", "theta4", "structure")
## 
## ####
## # Loop for Gumbel
## # The loop terminates, if 1000 structures are correctly identified
## 
## while(sum(results5_G[,"structure"]) < 1000){
##     data_hac5_G = rHAC(n, hac5_G)
## 
##     copula = estimate.copula(data_hac5_G, type = HAC_GUMBEL, method = ML)
## 
##     # Converts the structure to a string without dependency parameters
##     struc = tree2str(copula, theta=FALSE)
## 
##     # Check whether the structure is correct
##     correct = if((struc==struc1) | (struc==struc2) | (struc==struc3) |
##     (struc==struc4) | (struc==struc5) | (struc==struc6) | (struc==struc7) | (struc==struc8) |
##     (struc==struc9) | (struc==struc10) | (struc==struc11) | (struc==struc12) | (struc==struc13) |
##     (struc==struc14) | (struc==struc15) | (struc==struc16)){TRUE}else{FALSE}
## 
##     if(correct==TRUE){
##         results5_G = rbind(results5_G, c(get.params(copula, sort = TRUE, decreasing = TRUE), 1))
##     }else{
##         results5_G = rbind(results5_G, c(rep(NA,4), 0))
##     }
## }
## 
## ####
## # Loop for Clayton
## # The loop terminates, if 1000 structures are correctly identified
## 
## while(sum(results5_C[,"structure"]) < 1000){
##     data_hac5_C = rHAC(n, hac5_C)
## 
##     copula = estimate.copula(data_hac5_C, type = HAC_CLAYTON, method = ML)
## 
##     # Converts the structure to a string without dependency parameters
##     struc = tree2str(copula, theta=FALSE)
## 
##     # Check whether the structure is correct
##     correct = if((struc==struc1) | (struc==struc2) | (struc==struc3) |
##     (struc==struc4) | (struc==struc5) | (struc==struc6) | (struc==struc7) | (struc==struc8) |
##     (struc==struc9) | (struc==struc10) | (struc==struc11) | (struc==struc12) | (struc==struc13) |
##     (struc==struc14) | (struc==struc15) | (struc==struc16)){TRUE}else{FALSE}
## 
##     if(correct==TRUE){
##         results5_C = rbind(results5_C, c(get.params(copula, sort = TRUE, decreasing = TRUE), 1))
##     }else{
##         results5_C = rbind(results5_C, c(rep(NA,4), 0))
##     }
## }
## 
## ####
## # Remove the initial row and read the results for the 5-dimensional fully nested models
## results5_G = results5_G[-1,]
## summary(results5_G)
## apply(results5_G, 2, sd)
## 
## results5_C = results5_C[-1,]
## summary(results5_C)
## apply(results5_C, 2, sd)
## 
## #####Table 5---------------------------------------------------------------------------------------------
## 
## library("HAC")
## 
## ####
## # Construct the 5-dimensional hac objects for Gumbel and Clayton
## 
## hac5_G = hac(type = HAC_GUMBEL, tree = list(list("X1", "X2", tau2theta(2/3)), list("X3", "X4", tau2theta(1/3)), "X5", tau2theta(1/9)))
## hac5_C = hac(type = HAC_CLAYTON, tree = list(list("X1", "X2", tau2theta(2/3, type = HAC_CLAYTON)), list("X3", "X4", tau2theta(1/3, type = HAC_CLAYTON)), "X5", tau2theta(1/9, type = HAC_CLAYTON)))
## 
## ####
## # Define all structures, which correspond to the same copula as hac5_G and hac5_C respectively and save them as string, i.e., struc1, struc2, ...
## 
## hac5_G1 = hac(tree = list(list("X1", "X2", 4), list("X3", "X4", 2.5),"X5",1.1))
## struc1 = tree2str(hac5_G1, theta = FALSE)
## hac5_G2 = hac(tree = list(list("X2", "X1", 4), list("X3", "X4", 2.5),"X5",1.1))
## struc2 = tree2str(hac5_G2, theta = FALSE)
## hac5_G3 = hac(tree = list(list("X1", "X2", 4), list("X4", "X3", 2.5),"X5",1.1))
## struc3 = tree2str(hac5_G3, theta = FALSE)
## hac5_G4 = hac(tree = list(list("X2", "X1", 4), list("X4", "X3", 2.5),"X5",1.1))
## struc4 = tree2str(hac5_G4, theta = FALSE)
## hac5_G5 = hac(tree = list(list("X1", "X2", 4), "X5", list("X3", "X4", 2.5),1.1))
## struc5 = tree2str(hac5_G5, theta = FALSE)
## hac5_G6 = hac(tree = list(list("X2", "X1", 4), "X5", list("X3", "X4", 2.5),1.1))
## struc6 = tree2str(hac5_G6, theta = FALSE)
## hac5_G7 = hac(tree = list(list("X1", "X2", 4), "X5", list("X4", "X3", 2.5),1.1))
## struc7 = tree2str(hac5_G7, theta = FALSE)
## hac5_G8 = hac(tree = list(list("X2", "X1", 4), "X5", list("X4", "X3", 2.5),1.1))
## struc8 = tree2str(hac5_G8, theta = FALSE)
## hac5_G9 = hac(tree = list("X5", list("X1", "X2", 4), list("X3", "X4", 2.5),1.1))
## struc9 = tree2str(hac5_G9, theta = FALSE)
## hac5_G10 = hac(tree = list("X5", list("X2", "X1", 4), list("X3", "X4", 2.5),1.1))
## struc10 = tree2str(hac5_G10, theta = FALSE)
## hac5_G11 = hac(tree = list("X5", list("X1", "X2", 4), list("X4", "X3", 2.5),1.1))
## struc11 = tree2str(hac5_G11, theta = FALSE)
## hac5_G12 = hac(tree = list("X5", list("X2", "X1", 4), list("X4", "X3", 2.5),1.1))
## struc12 = tree2str(hac5_G12, theta = FALSE)
## hac5_G13 = hac(tree = list(list("X3", "X4", 2.5), list("X1", "X2", 4), "X5",1.1))
## struc13 = tree2str(hac5_G13, theta = FALSE)
## hac5_G14 = hac(tree = list(list("X3", "X4", 2.5), list("X2", "X1", 4), "X5",1.1))
## struc14 = tree2str(hac5_G14, theta = FALSE)
## hac5_G15 = hac(tree = list(list("X3", "X4", 2.5), list("X1", "X2", 4), "X5",1.1))
## struc15 = tree2str(hac5_G15, theta = FALSE)
## hac5_G16 = hac(tree = list(list("X3", "X4", 2.5), list("X2", "X1", 4), "X5",1.1))
## struc16 = tree2str(hac5_G16, theta = FALSE)
## hac5_G17 = hac(tree = list("X5", list("X3", "X4", 2.5), list("X1", "X2", 4),1.1))
## struc17 = tree2str(hac5_G17, theta = FALSE)
## hac5_G18 = hac(tree = list("X5", list("X3", "X4", 2.5), list("X1", "X2", 4),1.1))
## struc18 = tree2str(hac5_G18, theta = FALSE)
## hac5_G19 = hac(tree = list("X5", list("X4", "X3", 2.5), list("X1", "X2", 4),1.1))
## struc19 = tree2str(hac5_G19, theta = FALSE)
## hac5_G20 = hac(tree = list("X5", list("X4", "X3", 2.5), list("X1", "X2", 4),1.1))
## struc20 = tree2str(hac5_G20, theta = FALSE)
## hac5_G21 = hac(tree = list(list("X3", "X4", 2.5), "X5", list("X1", "X2", 4),1.1))
## struc21 = tree2str(hac5_G21, theta = FALSE)
## hac5_G22 = hac(tree = list(list("X3", "X4", 2.5), "X5", list("X1", "X2", 4),1.1))
## struc22 = tree2str(hac5_G22, theta = FALSE)
## hac5_G23 = hac(tree = list(list("X4", "X3", 2.5), "X5", list("X1", "X2", 4),1.1))
## struc23 = tree2str(hac5_G23, theta = FALSE)
## hac5_G24 = hac(tree = list(list("X4", "X3", 2.5), "X5", list("X1", "X2", 4),1.1))
## struc24 = tree2str(hac5_G24, theta = FALSE)
## 
## ####
## # Define the initial matrices for the results
## results5_G = matrix(0, nrow = 1, ncol = 4)
## struc5_G = matrix(0, nrow = 1, ncol = 1)
## results5_G_FML = matrix(0, nrow = 1, ncol = 3)
## 
## results5_C = matrix(0, nrow = 1, ncol = 4)
## struc5_C = matrix(0, nrow = 1, ncol = 1)
## results5_C_FML = matrix(0, nrow = 1, ncol = 3)
## 
## colnames(results5_G) = c("theta1", "theta2", "theta3", "structure")
## colnames(struc5_G) = c("realizedStruc")
## colnames(results5_G_FML) = c("theta1", "theta2", "theta3")
## 
## colnames(results5_C) = c("theta1", "theta2", "theta3", "structure")
## colnames(struc5_C) = c("realizedStruc")
## colnames(results5_C_FML) = c("theta1", "theta2", "theta3")
## 
## ####
## # Loop for Gumbel
## # The loop terminates, if 1000 structures are correctly identified
## 
## while(sum(results5_G[,"structure"]) < 1000){
##     data_hac5_G = rHAC(n, hac5_G)
## 
##     copula = estimate.copula(data_hac5_G, type = HAC_GUMBEL, method = RML, epsilon = 0.15)
## 
##     # Converts the structure to a string without dependency parameters
##     struc = tree2str(copula, theta=FALSE)
## 
##     # Check whether the structure is correct
##     correct = if((struc==struc1) | (struc==struc2) | (struc==struc3) |
##     (struc==struc4) | (struc==struc5) | (struc==struc6) | (struc==struc7) | (struc==struc8) |
##     (struc==struc9) | (struc==struc10) | (struc==struc11) | (struc==struc12) | (struc==struc13) |
##     (struc==struc14) | (struc==struc15) | (struc==struc16) | (struc==struc17) | (struc==struc18) |
##     (struc==struc19) | (struc==struc20) | (struc==struc21) | (struc==struc22) | (struc==struc23) |
##     (struc==struc24)){TRUE}else{FALSE}
## 
##     if(correct == TRUE){
##         results5_G = rbind(results5_G, matrix(c(get.params(copula, sort = TRUE, decreasing = TRUE), 1), nrow = 1, ncol =
##         4))
## 
##         # Reestimated the model with full ML
##         copula_FML = estimate.copula(data_hac5_G, type = HAC_GUMBEL, method = FML, hac = copula)
##         results5_G_FML = rbind(results5_G_FML, get.params(copula_FML, sort = TRUE, decreasing = TRUE))
##     }else{
##         results5_G = rbind(results5_G, matrix(c(rep(NA,3), 0), nrow = 1, ncol = 4))
##     }
##     struc5_G=rbind(struc5_G, struc)
## }
## 
## ####
## # Loop for Clayton
## # The loop terminates, if 1000 structures are correctly identified
## 
## while(sum(results5_C[,"structure"]) < 1000){
##     data_hac5_C = rHAC(n, hac5_C)
## 
##     copula = estimate.copula(data_hac5_C, type = HAC_CLAYTON, method = RML, epsilon = 0.2)
## 
##     # Converts the structure to a string without dependency parameters
##     struc = tree2str(copula, theta=FALSE)
## 
##     # Check whether the structure is correct
##     correct = if((struc==struc1) | (struc==struc2) | (struc==struc3) |
##     (struc==struc4) | (struc==struc5) | (struc==struc6) | (struc==struc7) | (struc==struc8) |
##     (struc==struc9) | (struc==struc10) | (struc==struc11) | (struc==struc12) | (struc==struc13) |
##     (struc==struc14) | (struc==struc15) | (struc==struc16) | (struc==struc17) | (struc==struc18) |
##     (struc==struc19) | (struc==struc20) | (struc==struc21) | (struc==struc22) | (struc==struc23) |
##     (struc==struc24)){TRUE}else{FALSE}
## 
##     if(correct == TRUE){
##         results5_C = rbind(results5_C, c(get.params(copula, sort = TRUE, decreasing = TRUE), 1))
## 
##         # Reestimated the model with full ML
##         copula_FML = estimate.copula(data_hac5_C, type = HAC_CLAYTON, method = FML, hac = copula)
##         results5_C_FML = rbind(results5_C_FML, get.params(copula_FML, sort = TRUE, decreasing = TRUE))
##     }else{
##         results5_C = rbind(results5_C, matrix(c(rep(NA,3), 0), nrow = 1, ncol = 4))
##         }
##     struc5_C=rbind(struc5_C, struc)
## }
## 
## ####
## # Remove the initial row and print the results for the 5-dimensional fully nested models
## results5_G = results5_G[-1,]
## struc5_G = struc5_G[-1,]
## results5_G_FML = results5_G_FML[-1,]
## results5_C = results5_C[-1,]
## struc5_C = struc5_C[-1,]
## results5_C_FML = results5_C_FML[-1,]
## 
## summary(na.omit(results5_G))
## apply(na.omit(results5_G), 2, sd)
## 
## summary(results5_G_FML)
## apply(results5_G_FML, 2, sd)
## 
## summary(results5_C)
## apply(na.omit(results5_C), 2, sd)
## 
## summary(results5_C_FML)
## apply(results5_C_FML, 2, sd)
## 
## ####
## # Compute the percentages of different classified structures
## # Our three cases for HAC_GUMBEL:
## # 1. True structures
## mean(struc5_G=="((X3.X4).(X1.X2).X5)") + mean(struc5_G=="((X1.X2).(X3.X4).X5)")
## 
## # Wrong aggregated models
## 
## round(mean(struc5_G=="((X1.X2).(X5.X3.X4))"),4)
## round(mean(struc5_G=="((X1.X2).((X3.X4).X5))"),4)
## round(mean(struc5_G=="((X1.X2).X5.X3.X4)"),4)
## 
## ####
## # Compute the percentages of different classified structures
## # Our three cases for HAC_CLAYTON:
## # 1. True structures
## 
## mean(struc5_C=="((X3.X4).(X1.X2).X5)") + mean(struc5_C=="((X1.X2).(X3.X4).X5)")
## 
## # Wrong aggregated models
## 
## mean(struc5_C=="(((X3.X4).(X1.X2)).X5)")
## mean(struc5_C== "((X1.X2).((X3.X4).X5))")
## 
## save(results5_G, struc5_G, results5_G_FML, file = "Sim_5_Gumbel_partially.Rdata")
## save(results5_C, struc5_C, results5_C_FML, file ="Sim_5_Clayton_partially.Rdata")


