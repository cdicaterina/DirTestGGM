### Scenario 1
### q = 6, n = 10 (as in Dawid & Lauritzen, 1993)
source("simulation_functions.R")
q <- 6
ybar <- 1:6

# null model - four edges
G0 <- UG(~ q1*q2 + q2*q3 + q4*q5 + q5*q6)

Lambda0 <- G0*0.5
diag(Lambda0) <- 1

Sigma0 <- solve(Lambda0)

rownames(Sigma0) <- colnames(Sigma0) <- paste("q", 1:q, sep = '')
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')

# unsaturated larger chordal model - seven edges
G1 <- UG(~ q1*q2 + q2*q3 + q2*q4 + q2*q5 + q3*q5 + q4*q5 + q5*q6)
rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')

p <- sum(G1)/2 + q
p0 <- sum(G0)/2 + q
d <- p - p0

nsim <- 100000
ncores <- 20

seed <- 123

n <- 8
wid.n <- 15

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

res_small <- simu_par(data.sim, wid = wid.n, cores = ncores)
save(res_small, file = "dir_small_res.rda")

