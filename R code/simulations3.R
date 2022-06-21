### Scenario 3
### Block diagonal concentration matrix
source("simulation_functions.R")
# concentration matrix Lambda (Omega in the paper) with q = 50
q <- 50
ybar <- rep(1:(q/5), 5)

# null model - two blocks (with 25 + 25 variables)
Sigma01 <- matrix(0.5, 25, 25)
Sigma02 <- matrix(0.5, 25, 25)

diag(Sigma01) <- rep(1, 25)
diag(Sigma02) <- rep(1, 25)

Sigma0 <- bdiag(Sigma01, Sigma02)
Sigma0 <- as.matrix(Sigma0)
Lambda0 <- solve(Sigma0)
Lambda0 <- (Lambda0 + t(Lambda0))/2
rownames(Sigma0) <- colnames(Sigma0) <- paste("q", 1:q, sep = '')
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')

G0 <- Lambda0
G0[Lambda0 == 0] <- 0
G0[Lambda0 != 0] <- 1
diag(G0) <- 0
rownames(G0) <- colnames(G0) <- paste("q", 1:q, sep = '')

# unsaturated model under H1
G1 <- G0
G1[16:25, 26:50] <- 1
G1[26:50, 16:25] <- 1
diag(G1) <- 0
rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')

p <- sum(G1)/2 + q
p0 <- sum(G0)/2 + q
d <- p - p0

nsim <- 100000
ncores <- 20

# maximal clique size: 35 < n

for(n in c(40, 60, 90, 120)) {
  if(n == 40) wid.n <- 50
  if(n == 60) wid.n <- 15
  if(n == 90) wid.n <- 10
  if(n == 120) wid.n <- 8
  
  data.sim <- as.list(numeric(nsim))
  set.seed(123)
  
  # generate Gaussian samples
  for(i in 1:nsim) {
    y <- rmvnorm(n, mean = ybar, sigma = Sigma0)
    colnames(y) <- colnames(Sigma0)
    data.sim[[i]] <- list(id = i, y = y)
  }
  
  attr(data.sim, "n") <- n
  attr(data.sim, "d") <- d
  attr(data.sim, "q") <- q
  attr(data.sim, "mu0") <- ybar
  attr(data.sim, "Sigma0") <- Sigma0
  attr(data.sim, "G0") <- G0
  attr(data.sim, "G1") <- G1
  
  res <- simu_par(data.sim, wid = wid.n, cores = ncores)
  save(res, file = paste0("dir_block_res_n", n, ".rda"))
}

