### Scenario 2
### Markovian dependence (as in Davison et al., 2014)
source("simulation_functions.R")

nsim <- 100000
ncores <- 20

n <- 60
# set q
q <- 11
wid.q <- 15

# load cow data
cow <- read.table("cowdata.txt", skip = 7, header = TRUE)
names(cow)[-(1:2)] <- paste("q", 1:11, sep = '')
cow.times <- cow[, -(1:2)]
ybar <- apply(cow.times, 2, mean)

# null model: MD(1)
Lambda0 <- matrix(0, nrow = q, ncol = q)
Lambda0[seq(2, q^2, by = q + 1)] <- rep(0.5, q - 1)
Lambda0 <- Lambda0 + t(Lambda0)
diag(Lambda0) <- rep(1, q)
Lambda0 <- diag(1:q) %*% Lambda0 %*% diag(1:q)
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')

G0 <- Lambda0
G0[Lambda0 == 0] <- 0
G0[Lambda0 != 0] <- 1
diag(G0) <- 0

model0 <- fitConGraph(G0, S = cov(cow.times), n = n)
Sigma0 <- model0$Shat

p0 <- sum(G0)/2 + q

# unsaturated alternative model: MD(m)
# set m < q - 1:
for(m in c(2, 3, 6, 9)) {
  G1 <- matrix(1, nrow = q, ncol = q)
  for(i in 1:(q - m - 1)) G1[i, (m + 1 + i):q] <- G1[(m + 1 + i):q, i] <- diag(G1) <- 0
  rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')
  p <- sum(G1)/2 + q
  d <- p - p0
  
  data.sim <- as.list(numeric(nsim))
  set.seed(123)
  
  # generate Gaussian sample
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
  
  res <- simu_par(data.sim, wid = wid.q, cores = ncores)
  save(res, file = paste0("dir_q11_res_MD", m, ".rda"))
}

rm(data.sim)
rm(res)

# set q
q <- 30
wid.q <- 15
ybar <- rep(1:(q/3), 3)

# null model: MD(1)
Lambda0 <- matrix(0, nrow = q, ncol = q)
Lambda0[seq(2, q^2, by = q + 1)] <- rep(0.5, q - 1)
Lambda0 <- Lambda0 + t(Lambda0)
diag(Lambda0) <- rep(1, q)
Lambda0 <- diag(rep(1:(q/10), 10)) %*% Lambda0 %*% diag(rep(1:(q/10), 10))
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')
Sigma0 <- solve(Lambda0)

G0 <- Lambda0
G0[Lambda0 == 0] <- 0
G0[Lambda0 != 0] <- 1
diag(G0) <- 0

p0 <- sum(G0)/2 + q

# unsaturated alternative model: MD(m)
# set m < q - 1:
for(m in c(2, 9, 18, 28)) {
  G1 <- matrix(1, nrow = q, ncol = q)
  for(i in 1:(q - m - 1)) G1[i, (m + 1 + i):q] <- G1[(m + 1 + i):q, i] <- diag(G1) <- 0
  rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')
  p <- sum(G1)/2 + q
  d <- p - p0
  
  data.sim <- as.list(numeric(nsim))
  set.seed(123)
  
  # generate Gaussian sample
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
  
  res <- simu_par(data.sim, wid = wid.q, cores = ncores)
  save(res, file = paste0("dir_q30_res_MD", m, ".rda"))
}

rm(data.sim)
rm(res)

# set q
q <- 50
wid.q <- 25
ybar <- rep(1:(q/5), 5)

# null model: MD(1)
Lambda0 <- matrix(0, nrow = q, ncol = q)
Lambda0[seq(2, q^2, by = q + 1)] <- rep(0.5, q - 1)
Lambda0 <- Lambda0 + t(Lambda0)
diag(Lambda0) <- rep(1, q)
Lambda0 <- diag(rep(1:(q/10), 10)) %*% Lambda0 %*% diag(rep(1:(q/10), 10))
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')
Sigma0 <- solve(Lambda0)

G0 <- Lambda0
G0[Lambda0 == 0] <- 0
G0[Lambda0 != 0] <- 1
diag(G0) <- 0

p0 <- sum(G0)/2 + q

# unsaturated alternative model: MD(m)
# set m < q - 1:
for(m in c(2, 16, 32, 48)) {
  G1 <- matrix(1, nrow = q, ncol = q)
  for(i in 1:(q - m - 1)) G1[i, (m + 1 + i):q] <- G1[(m + 1 + i):q, i] <- diag(G1) <- 0
  rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')
  p <- sum(G1)/2 + q
  d <- p - p0
  
  data.sim <- as.list(numeric(nsim))
  set.seed(123)
  
  # generate Gaussian sample
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
  
  res <- simu_par(data.sim, wid = wid.q, cores = ncores)
  save(res, file = paste0("dir_q50_res_MD", m, ".rda"))
}
