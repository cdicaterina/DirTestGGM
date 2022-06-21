##################################   Simulations   ###################################
### Packages
library(matrixcalc)
library(ggm)
library(Matrix)
library(mvtnorm)
library(gRbase)
library(plyr)
library(doMC)

### Find tmin e tmax
tmin <- function(that, Sig, c = 15, max.iter = 100) {
  iter <- 0
  min.t <- that - c * Sig
  if(min.t < 0) min.t <- 0
  while((min.t >= 1) & (iter <= max.iter)) {
    iter <- iter + 1
    c <- c + 2
    min.t <- that - c * Sig
    if(min.t < 0) min.t <- 0
  }
  if(iter > max.iter) min.t <- 1 - 1e-06
  return(min.t)
}

tmax <- function(that, Sig, tsup, c = 15, max.iter = 100) {
  iter <- 0
  max.t <- that + c * Sig
  if(max.t > tsup) max.t <- tsup
  while((max.t <= 1) & (iter <= max.iter)) {
    iter <- iter + 1
    c <- c + 2
    max.t <- that + c * Sig
    if(max.t > tsup) max.t <- tsup
  }
  if(iter > max.iter) max.t <- 1 + 1e-06
  max.t
}

### Tests
# w (from marginal likelihood based on Sigmahat)
lrt <- function(model0, model1, n) {
  if (!(is.list(model1))) {
    warning("model 1 taken to be saturated")
    app <- model1    
    model1 <- list()
    model1$Shat <- app
    model1$df <- 0
  }
  Sigma0 <- model0$Shat
  Sigma1 <- model1$Shat
  d <- model0$df - model1$df
  # statistic
  W <- (n - 1) * (determinant(Sigma0, logarithm = TRUE)$modulus[1] - 
                    determinant(Sigma1, logarithm = TRUE)$modulus[1])
  pvalue <- pchisq(W, d, lower.tail = FALSE)
  list(Woss = W, pvalue = pvalue)
}

# vech - uppertriangular part of symmetric matrix A
vech_ord <- function(A) {
  if(nrow(A) !=  ncol(A)) stop("matrix not square")
  if(!is.symmetric.matrix(A)) warning("matrix not symmetric")
  
  p <- nrow(A)
  pstar <- p * (p + 1)/2
  a <- matrix(NA, nrow = pstar, ncol = 1)
  index <- matrix(NA, nrow = pstar, ncol = 2)
  colnames(index) <- c("row", "col")
  var_names <- matrix(NA, nrow = pstar, ncol = 2)
  colnames(var_names) <- c("var1", "var2")
  
  # diagonal entries first
  for(i in 1:p) {
    a[i, 1] <- A[i, i]
    index[i, ] <- c(i, i)
    var_names[i, ] <- rep(colnames(A)[i], 2)
  }
  # upper-diagonal entries
  row_ind <- 1
  col_ind <- 2
  for(i in 1:(p * (p - 1)/2)) {
    a[i + p] <- A[row_ind, col_ind]
    index[i + p, ] <- c(row_ind, col_ind)
    var_names[i + p, ] <- c(colnames(A)[row_ind], colnames(A)[col_ind])
    if(col_ind ==  p) {
      row_ind <- row_ind + 1
      col_ind <- row_ind + 1
    }
    else 
      col_ind <- col_ind + 1
  }
  list(a = a, index = index, var_names = var_names)
}

# create vech(Lambda) = (lambda_k, lambda_h) = (lambda_k, 0)
lambda <- function(Lam) {
  vech <- vech_ord(Lam)
  lam <- vech$a
  lam[abs(lam) < 1e-6] <- 0    # larger threshold: 1e-6
  k <- lam !=  0
  v <- sum(k)
  lambda <- c(lam[k], lam[!k])
  index <- rbind(vech$index[k, ], vech$index[!k, ])
  var_names <- rbind(vech$var_names[k, ], vech$var_names[!k, ])
  list(lambda = lambda, index = index, var_names = var_names, v = v)
}

# Isserlis matrix
iss_kk <- function(A, index_k, v) {
  if(nrow(A) !=  ncol(A)) stop("matrix not square")
  if(!is.symmetric.matrix(A)) warning("matrix not symmetric")
  
  row <- index_k[1:v, 1]
  col <- index_k[1:v, 2]
  
  iss_kk <- A[row, row] * A[col, col] + A[row, col] * A[col, row]
  name <- paste0(row, ", ", col)
  dimnames(iss_kk) <- list(name, name)
  return(iss_kk)
}

# log determinant of Iss(Sigmamat)_kk with non-chordal option
logdet_iss_kk_f <- function(Sigmamat, GM, k_ind, length_k) {
  p <- nrow(Sigmamat)
  adj1 <- rip(GM)
  if(is.null(adj1)) {
    # non-chordal graph: std computation of Isserlis
    Iss_kk <- iss_kk(Sigmamat, index_k = k_ind, v = length_k)
    ldI <- tryCatch(determinant(Iss_kk)$modulus, error = function(e) 
      print("error ld_Iss_std"))
  }
  else {
    # chordal graph
    cliq <- adj1$cliques
    sep <- adj1$separators
    k <- length(cliq)
    s <- length(sep) - 1
    
    dets_c <- sapply(cliq, function(x) determinant(as.matrix(Sigmamat[x, x]))$modulus, 
                     simplify = TRUE)
    dets_s <- sapply(sep, function(x) determinant(as.matrix(Sigmamat[x, x]))$modulus, 
                     simplify = TRUE)[-1]
    
    card_c <- sapply(cliq, length)
    card_s <- sapply(sep[-1], length)
    
    ldI <- tryCatch(p * log(2) + sum((card_c + rep(1, k)) * dets_c) - 
                      sum((card_s + rep(1, s)) * dets_s), error = function(e) 
                        print("error ld_Iss"))
  }
  return(ldI)
}

# unsaturated model1 vs nested model0: Directional p-value
dir_p <- function(model0, model1, n, wid = 15, G1, plot = FALSE, npts = 100) {
  g.monot <- FALSE
  if (!(is.list(model1))) {
    app <- model1    
    model1 <- list()
    model1$Shat <- app
    model1$df <- 0
  }
  Sigma1 <- model1$Shat
  Sigma0 <- model0$Shat
  d <- model0$df - model1$df
    
  Lambda1 <- solve(Sigma1)
  if(!is.symmetric.matrix(Lambda1)) Lambda1 <- (Lambda1 + t(Lambda1))/2
    
  # t_sup analytic
  eigenvals <- eigen(solve(Sigma0) %*% Sigma1, only.values = TRUE)$values
  eigenvals <- as.double(eigenvals)
  t_sup <- 1/(1 - min(eigenvals))
    
  # 2nd term gbar
  gbar2 <- function(t) {
      (n - 1)/2 * sum(log(1 - t + t * eigenvals))
  }
    
  adj1 <- rip(G1)
  if(is.null(adj1)) {
    lambda1 <- lambda(Lambda1)
    Iss0_kk <- iss_kk(Sigma0, lambda1$index, lambda1$v)
    detIss0_kk <- determinant(Iss0_kk)$modulus
    # non-chordal graph: std computation of Isserlis matrix
    logdet_iss_kk <- function(t, k_ind, length_k) {
      Sigma <- Sigma0 + t * (Sigma1 - Sigma0)
      Iss_kk <- iss_kk(Sigma, index_k = k_ind, v = length_k)
      ldI <- tryCatch(determinant(Iss_kk)$modulus - detIss0_kk, 
                        error = function(e) print("error ld_Iss_std"))
      return(ldI)
    }
    gbar <- function(t) {
      ldI <- logdet_iss_kk(t, lambda1$index, lambda1$v)
      if(is.character(ldI)) return(NA)
      else return((d - 1) * log(t) + gbar2(t) - 0.5 * ldI)
    }
  }
  # chordal graph
  else {
    cliq1 <- adj1$cliques
    sep1 <- adj1$separators
    k1 <- length(cliq1)
    s1 <- length(sep1) - 1
      
    card1_c <- sapply(cliq1, length)
    card1_s <- sapply(sep1[-1], length)
      
    dets_c0 <- sapply(cliq1, function(x) determinant(as.matrix(Sigma0[x, x]))$modulus, 
                        simplify = TRUE)
    dets_s0 <- sapply(sep1, function(x) determinant(as.matrix(Sigma0[x, x]))$modulus, 
                        simplify = TRUE)[-1]
    # log determinant of Iss(Sigmahat(t))_kk
    logdet_iss_kk <- function(t, cliq, sep, k, s, dets0_c, dets0_s, card_c,
                                card_s) {
      Sigma <- Sigma0 + t * (Sigma1 - Sigma0)
      dets_c <- sapply(cliq, function(x) determinant(as.matrix(Sigma[x, x]))$modulus, 
                         simplify = TRUE)
      dets_s <- sapply(sep, function(x) determinant(as.matrix(Sigma[x, x]))$modulus, 
                         simplify = TRUE)[-1]
        
      ldI <- tryCatch(sum((card_c + rep(1, k)) * dets_c) - sum((card_s + 
                       rep(1, s)) * dets_s) - sum((card_c + rep(1, k)) * dets0_c) + 
                          sum((card_s + rep(1, s)) * dets0_s), 
                        error = function(e) print("error ld_Iss"))
      return(ldI)
    }
    gbar <- function(t) {
      ldI <- logdet_iss_kk(t, cliq1, sep1, k1, s1, dets_c0, dets_s0,
                           card1_c, card1_s)
      if(is.character(ldI)) return(NA)
      else return((d - 1) * log(t) + gbar2(t) - 0.5 * ldI)
    }
  }
    
  t_hat <- nlminb(1, function(x) -gbar(x), lower = 0, upper = t_sup)$par
  if(t_hat == t_sup) {
    g.monot <- TRUE
    gbarhat <- 0
    t_min <- 0
    t_max <- t_sup
  }
  else {
    gbarhat <- gbar(t_hat)
    # analytical -2nd derivative (incomplete)
    information <- function(t){
      (d - 1)/t^2 + 0.5 * (n - 1) * sum((1 - eigenvals)^2/(1 - t + t * eigenvals)^2) # + ...
    }
    sig <- 1/sqrt(information(t_hat))
      
    t_min <- tmin(t_hat, sig, c = wid)
    t_max <- tmax(t_hat, sig, t_sup, c = wid)
  }
    
  ff <- function(t) exp(gbar(t) - gbarhat)
  ff.v <- Vectorize(ff)
    
  up <- tryCatch(integrate(ff.v, lower = 1, upper = t_max),
                   error = function(e) print("error int sup"))
  down <- tryCatch(integrate(ff.v, lower = t_min, upper = 1),
                     error = function(e) print("error int inf"))
  pval <- up$value/(down$value + up$value)
  list(pvalue = pval, g.monot = g.monot)
}

# unsaturated model1 vs nested model0: w* and w**
Wstar <- function(model0, model1, n, G1) { 
  if (!(is.list(model1))) {
    app <- model1    
    model1 <- list()
    model1$Shat <- app
    model1$df <- 0
  }
  
  Sigma1 <- model1$Shat
  Sigma0 <- model0$Shat
  d <- model0$df - model1$df
    
  Lambda1 <- solve(Sigma1)
  if(!is.symmetric.matrix(Lambda1)) Lambda1 <- (Lambda1 + t(Lambda1))/2
  Lambda0 <- solve(Sigma0)
  if(!is.symmetric.matrix(Lambda0)) Lambda0 <- (Lambda0 + t(Lambda0))/2
    
  # lambda = vech(Lambda) calculation
  lambda1 <- lambda(Lambda1)  
  lambda1k <- lambda1$lambda[1:lambda1$v]
  lambda0k <- Lambda0[lambda1$index[1:lambda1$v, ]]
  lambda0k[abs(lambda0k) < 1e-6] <- 0 
  # sigma = vech(Sigma) calculation
  sigma1k <- Sigma1[lambda1$index[1:lambda1$v, ]]
  sigma0k <- Sigma0[lambda1$index[1:lambda1$v, ]]
    
  p <- nrow(Sigma1)
  J_kk <- diag(c(rep(1, p), rep(2, lambda1$v - p)), nrow = p + lambda1$v - p,
                 ncol = p + lambda1$v - p)
    
  jpsi1 <- 0.25 * (J_kk %*% iss_kk(Sigma0, lambda1$index, lambda1$v) %*% J_kk)
    
  # LRT
  w <- lrt(model0, model1, n)
    
  # gamma
  num <- (d/2) * log(0.25 * (n - 1) * t(J_kk %*% (sigma0k - sigma1k)) %*% 
                         solve(jpsi1) %*% (J_kk %*% (sigma0k - sigma1k)))
  den <- log(0.5* (n - 1) * t(lambda1k - lambda0k) %*% (J_kk %*% 
                                                            (sigma0k - sigma1k)))
  ld_I0 <- logdet_iss_kk_f(Sigma0, G1, lambda1$index, lambda1$v)
  ld_I1 <- logdet_iss_kk_f(Sigma1, G1, lambda1$index, lambda1$v)
    
  if(is.character(ld_I0) | is.character(ld_I1)) wstar1 <- wstar2 <- 
      p1 <- p2 <- NA
  else {
    fact <- 0.5 * (ld_I0 - ld_I1)
    attr(fact, "logarithm") <- NULL
    log_gamma <- (num - den - ((d/2) - 1) * log(w$Woss) + fact)
    wstar1 <- w$Woss * (1 - log_gamma/w$W)^2
    wstar2 <- w$Woss - 2 * log_gamma
      
    p1 <- pchisq(wstar1, d, lower.tail = FALSE)
    p2 <- pchisq(wstar2, d, lower.tail = FALSE)
  }
  list(wstar1 = wstar1, p.value1 = p1, wstar2 = wstar2, p.value2 = p2)
}

### Simulations
covselect.simu1 <- function(data, G0, G1, wid = 15, trace = TRUE) {
  if(trace) if(data$id %% 10^3 == 0) cat(data$id, '\n')
  y <- as.data.frame(data$y)
  n <- nrow(y)
  S <- cov(y)
  
  # lrt 1st order
  model0s <- fitConGraph(G0, S, n)
  model1s <- fitConGraph(G1, S, n)
  
  lrt.p <- lrt(model0s, model1s, n = n)$pvalue
  
  # directional p-values
  dir <- dir_p(model0s, model1s, n = n, wid = wid, G1 = G1)
  dir.p <- dir$pvalue
  g.monot <- dir$g.monot
  
  # Skovgaard's p-values
  app <- Wstar(model0s, model1s, n = n, G1 = G1)
  sk1.p <- app$p.value1
  sk2.p <- app$p.value2
  
  data.frame(id = data$id, first.order = lrt.p, skovgaard1 = sk1.p,
             skovgaard2 = sk2.p, directional = dir.p, g.monot = g.monot, wid = wid)
}

simu_par <- function(data_all, cores = 1, wid, trace = TRUE, ncores = 2) {
  registerDoMC(ncores)
  
  Nsim <- length(data_all)
  n <- attr(data_all, "n")
  G0 <- attr(data_all, "G0")
  G1 <- attr(data_all, "G1")
  
  out <- ldply(.data = data_all, .fun = covselect.simu1, trace = trace, 
               .parallel = TRUE, G0 = G0, G1 = G1, wid = wid, .inform = TRUE, .id = id)
  out
}
