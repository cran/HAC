# estimate.r ############################################################################################################# 
# FUNCTION:               	DESCRIPTION: 
# .logLik_oneHierachy			  log-density of a HAC with one nested HAC (internal function)
# .pdf.1st.part             1st-part of the log-pdf (internal function)
# .b                        auxilary function for computation of density (internal function)
# .recB                     recurisve computation of bell-polynomials (internal function)
# .d.phi.inv.circ.phi       d-th derivative of .phi.inv.circ.phi (internal function)
# .d.phi.inv.circ.phi.aux   auxilary function for computation of d-th derivative of .phi.inv.circ.phi (internal function) 
# .pdf.2nd.part             2nd-part of the log-pdf (internal function)
##########################################################################################################################

.logLik_oneHierachy = function(U, tree, type){
    sum(.pdf.1st.part(U, tree, type) + .pdf.2nd.part(U, tree, type))  
}

#----------------------------------------------------------------------------------------

.pdf.1st.part = function(U, tree, type){
  d0 = length(unlist(tree[which(sapply(tree, is.character))]))
  rtheta = tree[[length(tree)]]
  
  ntree = which(sapply(tree, is.list))
  nvars = unlist(tree[[ntree]][which(sapply(tree[[ntree]], is.character))]); d1 = length(nvars)
  ntheta = tree[[ntree]][[length(tree[[ntree]])]]
  
  U0  = phi.inv(x = .cop.cdf(sample = U, tree = tree, type = type), theta = rtheta, type = type)
  U1 = rowSums(matrix(phi.inv(U[, nvars], theta = ntheta, type = type), ncol = d1))
  
  ind = 1:d1 + d0
  if(type == 1){
    phi = sapply(ind, copGumbel@absdPsi, t = U0, theta = rtheta)
  }
  if(type == 3){
    phi = sapply(ind, copClayton@absdPsi, t = U0, theta = rtheta)
  }
  if(type == 5){
    phi = sapply(ind, copFrank@absdPsi, t = U0, theta = rtheta)
  }
  if(type == 7){
    phi = sapply(ind, copJoe@absdPsi, t = U0, theta = rtheta)
  }
  if(type == 9){
    phi = sapply(ind, copAMH@absdPsi, t = U0, theta = rtheta)
  }
  log(rowSums(abs(phi * .b(U1 = U1, d1 = d1, rtheta = rtheta, ntheta = ntheta, type = type))))
}

#--------------------------------------------------------------------------------

.b = function(U1, d1, rtheta, ntheta, type){
  b = derivs = matrix(0, nrow = length(U1), ncol = d1)
  
  for(k in 1:d1){
      derivs[, k] = .d.phi.inv.circ.phi(x = U1, 
                                        rtheta = rtheta, 
                                        ntheta = ntheta, 
                                        k = k, 
                                        type = type)
  }
  
  J = list(); length(J) = d1
  for(i in 1:d1){
    J[[i]] = list(); length(J[[i]]) = d1 - 1
    for(j in 1:(d1 - 1)){
      J[[i]][[j]] = list(); length(J[[i]][[j]]) = 2
      J[[i]][[j]][[1]] = j:(i - 1) -> ind
      J[[i]][[j]][[2]] = choose(i, ind)
    }
  }
  for(j in 1:d1){
    b[, j] = apply(derivs, 1, FUN = .recB, n = d1, k = j - 1, J = J) / factorial(j)
  }
  b
}

#--------------------------------------------------------------------------------

.recB = function(n, k, x, J){
  if(k < 1){
    x[n]
  }else{
    sum(J[[n]][[k]][[2]] * x[n - J[[n]][[k]][[1]]] * sapply(J[[n]][[k]][[1]], .recB, k = k - 1, x = x, J = J))
  }
}

#--------------------------------------------------------------------------------

.d.phi.inv.circ.phi = function(x, rtheta, ntheta, k, type){
  alpha = rtheta / ntheta
  if(type == 1){
    return(prod(alpha - 1:k + 1) * x^(alpha - k))
  }
  if(type == 3){
    return(prod(alpha - 1:k + 1) * (1 + x)^(alpha - k))
  }
  #if(type == 5){
    #k = 0: -log(1 + exp(-rtheta) * (exp(rtheta - ( -log(1 - exp(-x) + exp(-ntheta - x))/ntheta) * rtheta) - 1)/(exp(-rtheta) - 1))
    #k = 0: -rtheta + log(-1 + exp(rtheta)) - log(1 - (1 - exp(-x) + exp(-x - ntheta))^(alpha))
  #  if(k == 1){
  #    -(alpha * (exp(ntheta) - 1) * (exp(-x - ntheta) - exp(-x) + 1)^alpha)/((-exp(ntheta) + exp(ntheta + x) + 1) * (1 - (exp(-x - ntheta) - exp(-x) + 1)^alpha))
  #  }else{
  #  a  
  #  }
  #}
  if(type == 7){
    #k = 0: -log(1 - (1 - exp(-x))^(alpha))
    if(k == 1){
      (alpha * (1 - exp(-x))^alpha)/((exp(x) - 1) * (1 - (1 - exp(-x))^alpha))
    }else{
      deriv = matrix(0, nrow = length(x), ncol = k)
      for(l in 0:(k - 1)){
        deriv[, l + 1] = choose(k - 1, l) * .d.phi.inv.circ.phi.aux(x, rtheta, ntheta, l, type)
      }
      return(rowSums(deriv)) 
    }
  }
  if(type == 9){
    if(k == 1){
      ((rtheta - 1) * exp(x))/(ntheta - rtheta + (rtheta - 1) * exp(x))
    }else{
      deriv = matrix(0, nrow = length(x), ncol = k)
      for(l in 0:(k - 1)){
        deriv[, l + 1] = choose(k - 1, l) * .d.phi.inv.circ.phi.aux(x, rtheta, ntheta, l, type)
      }
      return((rtheta - 1) * rowSums(exp(x) * deriv))
    }
  }
}

#--------------------------------------------------------------------------------

.d.phi.inv.circ.phi.aux = function(x, rtheta, ntheta, l, type){
  if(l > 0){
    J = list(); length(J) = l
    if(l > 1){
      for(i in 1:l){
        J[[i]] = list(); length(J[[i]]) = l - 1
        for(j in 1:(l - 1)){
          J[[i]][[j]] = list(); length(J[[i]][[j]]) = 2
          J[[i]][[j]][[1]] = j:(i - 1) -> ind
          J[[i]][[j]][[2]] = choose(i, ind)
        }
      }
    }
  }
  if(type == 5){
    alpha = rtheta / ntheta
    -(alpha * (exp(ntheta) - 1) * (exp(-x - ntheta) - exp(-x) + 1)^alpha)/((-exp(ntheta) + exp(ntheta + x) + 1) * (1 - (exp(-x - ntheta) - exp(-x) + 1)^alpha))
  }
  if(type == 7){
    alpha = rtheta / ntheta
    #1/((exp(x) - 1) * (1 - (1 - exp(-x))^alpha))
  }
  if(type == 9){
    if(l == 0){
      1 / (ntheta - rtheta + (rtheta - 1) * exp(x))
    }else{
      .in = .out = matrix(rep((rtheta - 1) * exp(x), l), nrow = length(x), ncol = l)
      for(ll in 1:l){
        .out[, ll] = (-1)^ll / (ntheta - rtheta + (rtheta - 1) * exp(x))^(ll + 1) * apply(.in, 1, FUN = .recB, n = l, k = ll - 1, J = J)
      }
      rowSums(.out)
    }
  }  
}

#--------------------------------------------------------------------------------

.pdf.2nd.part = function(U, tree, type){
  rvars = unlist(tree[which(sapply(tree, is.character))]); d0 = length(rvars)
  rtheta = tree[[length(tree)]]
  
  ntree = which(sapply(tree, is.list))
  nvars = unlist(tree[[ntree]][which(sapply(tree[[ntree]], is.character))]); d1 = length(nvars)
  ntheta = tree[[ntree]][[length(tree[[ntree]])]]
  if(type == 1){
    U0 = matrix(copGumbel@absdiPsi(U[, rvars], theta = rtheta, log = TRUE), ncol = d0)
    U1 = matrix(copGumbel@absdiPsi(U[, nvars], theta = ntheta, log = TRUE), ncol = d1)  
  }
  if(type == 3){
    U0 = matrix(copClayton@absdiPsi(U[, rvars], theta = rtheta, log = TRUE), ncol = d0)
    U1 = matrix(copClayton@absdiPsi(U[, nvars], theta = ntheta, log = TRUE), ncol = d1)
  }
  if(type == 5){
    U0 = matrix(copFrank@absdiPsi(U[, rvars], theta = rtheta, log = TRUE), ncol = d0)
    U1 = matrix(copFrank@absdiPsi(U[, nvars], theta = ntheta, log = TRUE), ncol = d1)
  }
  if(type == 7){
    U0 = matrix(copJoe@absdiPsi(U[, rvars], theta = rtheta, log = TRUE), ncol = d0)
    U1 = matrix(copJoe@absdiPsi(U[, nvars], theta = ntheta, log = TRUE), ncol = d1)
  }
  if(type == 9){
    U0 = matrix(copAMH@absdiPsi(U[, rvars], theta = rtheta, log = TRUE), ncol = d0)
    U1 = matrix(copAMH@absdiPsi(U[, nvars], theta = ntheta, log = TRUE), ncol = d1)
  }
  rowSums(cbind(U0, U1))
}