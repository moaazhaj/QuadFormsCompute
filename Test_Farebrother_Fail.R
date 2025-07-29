library(CompQuadForm)
library(bigQF)
source("spa.R")

ExampleParameters <- list(
  ### Example 1: Failure in Imhoff q > 35
  list(
    lambda  = exp(-(0:9)), 
    hr      = rep(1,   10),
    delta2  = rep(1,   10),
    q = seq(10, 50, length.out=5)
  ),
  list(
    lambda = c(1,2),
    hr     = c(1,1),
    delta2 = c(1,1e4),
    q = seq(15e3, 2e4, length.out=5)
  )
)

## Choose an example number
exNumber <-2

# Retrieve the chosen parameter set
params  <- ExampleParameters[[exNumber]]
lambda  <- params$lambda
delta2  <- params$delta2
hr      <- params$hr
q_vals   <- params$q


ccdf_farebrother1 <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-4, mode = 1)$Qq)

ccdf_farebrother2 <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-6, mode = 1)$Qq)

ccdf_farebrother3 <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-8, mode = 1)$Qq)

ccdf_farebrother4 <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-12, mode = 1)$Qq)

ccdf_davies <- sapply(q_vals, function(q) 
  davies(q, lambda, h = hr, delta = delta2, 
         sigma = 0, lim = 1e5, acc = 1e-5)$Qq)







###############################################################################
# Plotting code for the CCDF methods
###############################################################################
  
plot(q_vals, ccdf_farebrother1, type = "o", 
    # ylim = range(ccdf_farebrother1,ccdf_davies), 
     #ylim = c(1e-5,1e-1), 
     #xlim = c(0,120), 
     xlab = "q", 
     ylab = "CCDF", 
     lwd = 2, 
     col = "black",
     log = "y",
     main = paste("Comparison of CCDF  - Examples of Farebrother's Failure", exNumber))

lines(q_vals, ccdf_farebrother2,       col = "red",   lwd = 2, lty = 2)
lines(q_vals, ccdf_farebrother3,         col = "blue",  lwd = 1, lty = 3,pch = 20, type = "o")
lines(q_vals, ccdf_farebrother4,         col = "purple",lwd = 2, lty = 5)
lines(q_vals, ccdf_davies,  col = "green", lwd = 3, lty = 4)

# Add grid for better visualization on log scale
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")

# # Only plot Zhang, LPB, WF if they are not NA
# if(!all(is.na(ccdf_zhang))) {
#   lines(q_vals, ccdf_zhang,  col = "orange", lwd = 2, lty = 6)
# }

legend("bottomleft", 
       legend = c("Farebrother(acc = 1e-4)", "Farebrother(acc = 1e-6)", 
                  "Farebrother(acc = 1e-8)", "Farebrother(acc = 1e-12)", "Davies"),
       col = c("black", "red", "blue", 
               "purple", "green" ),
       lty = c(1,2,3,4,5),
       bty = "n", cex = 0.8)


