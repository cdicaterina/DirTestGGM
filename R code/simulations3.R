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

seed <- 123
for (n in c(60, 90, 120)) {
  res50p <- covselect.simu(mu = ybar, Sigmatrue = Sigma0, G0 = G0, 
                           G1 = G1, Nsim = nsim, n = n, seed = seed)
  save(res50p, file = paste("dir_block_res_p50n", n, ".rda", sep = ""))
}