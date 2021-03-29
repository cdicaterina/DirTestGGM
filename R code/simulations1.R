### Scenario 1
### q = 4, n = 7 (as in Eriksen, 1996)
source("simulation_functions.R")
q <- 4
ybar <- 1:4

# null model - two edges
G0 <- UG(~ q1*q2 + q3*q4)

# concentration matrix Lambda (Omega in the paper) with q = 4
Lambda0 <- G0*0.5
diag(Lambda0) <- 1

Sigma0 <- solve(Lambda0)
Sigma0 <- Sigma0/diag(Sigma0)
Lambda0 <- solve(Sigma0)

rownames(Sigma0) <- colnames(Sigma0) <- paste("q", 1:q, sep = '')
rownames(Lambda0) <- colnames(Lambda0) <- paste("q", 1:q, sep = '')

# unsaturated larger model - square
G1 <- UG(~ q1*q2 + q2*q3 + q3*q4 + q4*q1)
rownames(G1) <- colnames(G1) <- paste("q", 1:q, sep = '')

p <- sum(G1)/2 + q
p0 <- sum(G0)/2 + q
d <- p - p0

nsim <- 100000

seed <- 123

n <- 7

res_small <- covselect.simu(mu = ybar, Sigmatrue = Sigma0, G0 = G0, 
                            G1 = G1, Nsim = nsim, n = n, seed = seed)
save(res_small, file = paste("nondec_res_q4", n, ".rda", sep = ""))
