### Scenario 2
### Markovian dependence (as in Davison et al., 2014)
source("simulation_functions.R")
# set p (q in the paper), for instance:
p <- 50
n <- 60
ybar <- rep(1:(p/5), 5)

# null model MD(1)
# concentration matrix Lambda (Omega in the paper)
Lambda0 <- matrix(0, nrow = p, ncol = p)
Lambda0[seq(2, p^2, by = p + 1)] <- rep(0.5, p - 1)
Lambda0 <- Lambda0 + t(Lambda0)
diag(Lambda0) <- rep(1, p)
Lambda0 <- diag(rep(1:(p/10), 10)) %*% Lambda0 %*% diag(rep(1:(p/10), 10))
rownames(Lambda0) <- colnames(Lambda0) <- paste("p", 1:p, sep = '')
Sigma0 <- solve(Lambda0)
df <- length(which(Lambda0 == 0))/2

G0 <- Lambda0
G0[Lambda0 == 0] <- 0
G0[Lambda0 != 0] <- 1
diag(G0) <- 0

model0 <- fitConGraph(G0, S = Sigma0, n = n)

# unsaturated model MD(k)  (k is m in the paper)
# set k < p - 1, for instance:
k <- 32
G1 <- matrix(1, nrow = p, ncol = p)
for(i in 1:(p - k - 1)) G1[i, (k + 1 + i):p] <- G1[(k + 1 + i):p, i] <- diag(G1) <- 0
rownames(G1) <- colnames(G1) <- paste("p", 1:p, sep = '')

nsim <- 100000
Sigmatrue <- model0$Shat

seed <- 123
resMD32 <- covselect.simu(mu = ybar, Sigmatrue = Sigmatrue, G0 = G0, G1 = G1, 
                         Nsim = nsim, n = n, seed = seed)
save(resMD32, file = "dir_res_q50.rda")
