library(CompQuadForm)
library(bigQF)
source("spa.R")

ExampleParameters <- list(
  ### Example 1: Failure in Imhoff q > 35
  list(
    lambda  = c(1, 0.6^4),
    hr      = c(1,   1),
    delta2  = c(1.0, 7),
    q = seq(30, 40, length.out=50)
  )
)

## Choose an example number
exNumber <- 1

# Retrieve the chosen parameter set
params  <- ExampleParameters[[exNumber]]
lambda  <- params$lambda
delta2  <- params$delta2
hr      <- params$hr
q_vals   <- params$q


ccdf_imhof <- sapply(q_vals, function(q) 
  imhof(q, lambda, h = hr, delta = delta2, 
        epsabs = 1e-6, epsrel = 1e-8, limit = 1e6)$Qq)

ccdf_davies <- sapply(q_vals, function(q) 
  davies(q, lambda, h = hr, delta = delta2,  sigma = 0, lim = 5e7, acc = 1e-10)$Qq)

ccdf_spa <- sapply(q_vals, function(q) 1 - spa(q,lambda, hr, delta2))

ccdf_farebrother <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-12, mode = 1)$Qq)

ccdf_liu <- sapply(q_vals, function(q) liu(q, lambda, h = hr, delta = delta2))


###############################################################################
# Plotting code for the CCDF methods
###############################################################################
  
plot(q_vals, ccdf_imhof, type = "o", 
     ylim = c(1e-8,1e-4), 
     #xlim = c(0,120), 
     xlab = "q", 
     ylab = "CCDF", 
     lwd = 2, 
     col = "black",
     log = "y",
     main = paste("Comparison of CCDF  - Examples from Imhof", exNumber))

lines(q_vals, ccdf_davies,       col = "red",   lwd = 2, lty = 2)
lines(q_vals, ccdf_spa,          col = "blue",  lwd = 1, lty = 3,pch = 20, type = "o")
lines(q_vals, ccdf_farebrother,  col = "green", lwd = 3, lty = 4)
lines(q_vals, ccdf_liu,          col = "purple",lwd = 2, lty = 5)
# Add grid for better visualization on log scale
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")

# # Only plot Zhang, LPB, WF if they are not NA
# if(!all(is.na(ccdf_zhang))) {
#   lines(q_vals, ccdf_zhang,  col = "orange", lwd = 2, lty = 6)
# }

legend("topleft", 
       legend = c("Imhof", "Davies", "Saddlepoint", "Farebrother",
                  "Liu"),
       col = c("black", "red", "blue", "green", 
               "purple" ),
       lty = c(1,2,3,4,5),
       bty = "n", cex = 0.8)


