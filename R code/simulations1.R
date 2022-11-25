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

n <- 10
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

res_small_chord <- simu_par(data.sim, wid = wid.n, cores = ncores)
save(res_small_chord, file = "dir_small_res.rda")

# ##### Output
# ###################################    Output    ####################################
# ### Summary function
# summary.sim <- function(out, levels = c(0.01, 0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 
#                                         0.9, 0.95, 0.975, 0.99)) {
#   fo <- sapply(levels, function(x) mean(out$first.order < x))
#   sk1 <- sapply(levels, function(x) mean(out$skovgaard1 < x))
#   sk2 <- sapply(levels, function(x) mean(out$skovgaard2 < x))
#   dir2 <- sapply(levels, function(x) mean(out$directional < x))
#   
#   sim.se <- sqrt(levels * (1 - levels)/length(out$first.order))
#   res <- rbind(levels, fo, sk1, sk2, dir2, sim.se)
#   rownames(res) <- c("nominal", "first-order", "Skovgaard W*", 
#                      "Skovgaard W**", "directional", "sim-error")
#   res
# }
# 
# ### Figure 1
# library(igraph)
# 
# ug1 <- UG(~ A*B + B*C + B*D + B*E + C*E + D*E + E*F)
# graph1 <- graph.adjacency(ug1, mode = "undirected", diag = FALSE)
# 
# x <- tkplot(graph1, vertex.color = "lightgrey", vertex.label.color = "black", 
#             vertex.label.family = "Arial", vertex.size = 15, vertex.label.cex = 1.5)
# 
### Figure 2
# setwd("~/Dropbox/Directional/results")
setwd("~/")
load("dir_small_res.rda")

setwd("~/Dropbox/Directional/paper/ArXiv submission/figs")
pdf("fig0n8.pdf", width = 14, height = 5)
par(mfrow = c(1, 2), mai = c(1, 1, 0, 0.5))
old.pty <- par("pty")
par(pty = "s")
# 1st plot: empirical p-values
index <- seq(1, length(res_small_chord$first.order), by = 10)
x <- sort(res_small_chord$first.order)
x <- x[index]
plot(ppoints(x), x, type = "l",
     xlim= c(0, 0.1), lwd = 1.5, ylim= c(0, 0.1), xlab = "Uniform quantiles",
     ylab = "p-value", pty = "s", col = "chocolate2", lty = "dotdash", cex.lab = 1)
abline(0, 1, col ="grey")
x <- sort(res_small_chord$skovgaard1)
x <- x[index]
lines(ppoints(x), x, col = "seagreen3",
      lwd = 1.5, lty = "dashed")
x <- sort(res_small_chord$skovgaard2)
x <- x[index]
lines(ppoints(x), x, col = "seagreen4",
      lty = "longdash", lwd = 1.5)
x <- sort(res_small_chord$directional)
x <- x[index]
lines(ppoints(x), x, col = "cornflowerblue", lwd = 1.5)

par(pty = old.pty)

# 2nd plot: relative error for p-values lower than 0.10 (without first-order)
x <- sort(res_small_chord$directional)
x <- x[index]
plot(ppoints(x), (x - ppoints(x))/ppoints(x), xlim = c(0, 0.1), lwd = 1.5,
     ylab ="Relative error", xlab = "Nominal level", type ="l", ylim = c(-0.25, 0.15),
     cex.lab = 1, col = "cornflowerblue")
abline(h = 0, col = "grey")
x <- sort(res_small_chord$skovgaard1)
x <- x[index]
lines(ppoints(x), (x - ppoints(x))/ppoints(x),
      col = "seagreen3", lty = "dashed", lwd = 1.5)
x <- sort(res_small_chord$skovgaard2)
x <- x[index]
lines(ppoints(x), (x - ppoints(x))/ppoints(x), lwd = 1.5,
      col = "seagreen4", lty = "longdash")
par(mfrow = c(1, 1))
dev.off()

