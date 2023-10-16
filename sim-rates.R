## one simulation replication
## target parameter is alpha_1

epsilon <- 0.15 # bounding treatment probabilities away from zero
replication <- function(n, err) {
  
  h <- 2.2
  X1 <-  rbinom(n, 1, 0.5) # one confounder
  X2 <- runif(n) # another confounder
  
  tau.den <- plogis(-0.2 + 1.2 * X2 + 0.5 * X1)
  S <- rbinom(n, 1, tau.den) # 2 data sources (1 = trial and 0 = target)
  
  q.1 <- pmin(pmax(epsilon, plogis(0.3 - 0.8 * X2 + 0.9 * X1)),1 - epsilon)
  A <- rbinom(n, 1, ifelse(S == 1, q.1, 0))
  
  mu.0 <- 0.75 * (5.2 + X1 - 1.2 * X2)
  mu.11 <- 5.2 + 1.2 - 0.6 * X1 + X1 - 1.2 * X2
  mu.10 <- 5.2 + X1 - 1.2 * X2
  Y <- rnorm(n, ifelse(S == 0, mu.0, A * mu.11 + (1-A) * mu.10), 1)
  
  q.hat <- plogis(qlogis(q.1) + h * rnorm(n, 1/n^err, 1/n^err))
  tau.den.hat <- plogis(qlogis(tau.den) + h * rnorm(n, 1/n^err, 1/n^err))
  tau.hat <- (1 - tau.den.hat) / tau.den.hat
  mu.0.hat <- mu.0 + h * rnorm(n, 1/n^err, 1/n^err)
  mu.11.hat <- mu.11 + h * rnorm(n, 1/n^err, 1/n^err)
  mu.10.hat <- mu.10 + h * rnorm(n, 1/n^err, 1/n^err)
  
  r.A.hat <- mu.11.hat / mu.10.hat
  r.S.hat <- mu.0.hat / mu.10.hat
  
  ## plug-in estimator
  plug.l <- (1 - S) * r.A.hat * mu.0.hat / mean(1 - S)
  
  ## IF-based estimator
  IF.l <- ((1 - S) * r.A.hat * Y +
             S * tau.hat * r.S.hat *
             (A * (Y - mu.11.hat) / q.hat - 
                (1 - A) * r.A.hat * (Y - mu.10.hat) / (1 - q.hat))) / 
    mean(S == 0)
  
  return(c(plug = mean(plug.l), IF = mean(IF.l)))
  
}

set.seed(591)
R <- 5000
n <- 5000
err <- seq(0.10, 0.5, by = 0.05)

res <- sapply(1:length(err), function(l) {
  err.l <- err[l]
  sapply(1:R, FUN = function(j) {
    replication(n, err.l)
  })
}, simplify = "array")
res <- aperm(res, c(3, 2, 1))

### TRUE COUNTERFACTUAL MEAN CALCULATION
N <- 10000000
X1 <-  rbinom(n, 1, 0.5) # one confounder
X2 <- runif(n) # another confounder

S <- rbinom(n, 1, plogis(-0.2 + 1.2 * X2 + 0.5 * X1))
A <- rbinom(n, 1, ifelse(S == 1, pmin(pmax(epsilon, plogis(0.3 - 0.8 * X2 + 0.9 * X1)),1 - epsilon), 0))
Y <- rnorm(n, ifelse(S == 0, 0.75 * (5.2 + X1 - 1.2 * X2),
                     A * (5.2 + 1.2 - 0.6 * X1 + X1 - 1.2 * X2) + 
                       (1-A) * (5.2 + X1 - 1.2 * X2)), 1)
truth <- mean((5.2 + 1.2 - 0.6 * X1[S == 0] + X1[S == 0] - 1.2 * X2[S == 0]) / 
                (5.2 + X1[S == 0] - 1.2 * X2[S == 0]) * Y[S == 0])

rmse <- function(res) {
  sqrt(colMeans((res - truth)^2))
}

res.rmse <- cbind.data.frame(err, t(apply(res, 1, rmse)))

# par(mar = c(4,4,1.6,1), mfrow = c(3,1))
plot(res.rmse$err, res.rmse$plug, type = "b", xlim = c(0.10, 0.5), pch=20,
     ylim = c(0, max(c(res.rmse$IF, res.rmse$plug))), 
     lty = "dashed", col = "red", ylab = "RMSE", xlab = "Error Rate",
     main = expression(n ~ "=" ~ 5000))
lines(res.rmse$err, res.rmse$IF, type = "b", pch = 20)
legend("topright", lty = c("dashed","solid"), pch = 20,
       col = c("red", "black"), bty = "n",
       legend = c("Plug-in", "IF-based"))
