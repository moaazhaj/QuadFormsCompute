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
  ),
  list(
    lambda  = c(0.6, 0.3,0.1),
    hr      = c(1,  1 ,1),
    delta2  = c(0, 0, 0),
    q = seq(15, 17, length.out=50)
  )
)

## Choose an example number
exNumber <- 2

# Retrieve the chosen parameter set
params  <- ExampleParameters[[exNumber]]
lambda  <- params$lambda
delta2  <- params$delta2
hr      <- params$hr
q_vals   <- params$q


ccdf_imhof1 <- sapply(q_vals, function(q) 
  imhof(q, lambda, h = hr, delta = delta2, 
        epsabs = 1e-8, epsrel = 1e-8, limit = 1e6)$Qq)

ccdf_imhof2 <- sapply(q_vals, function(q) 
  imhof(q, lambda, h = hr, delta = delta2, 
        epsabs = 1e-9, epsrel = 1e-8, limit = 1e6)$Qq)

ccdf_imhof3 <- sapply(q_vals, function(q) 
  imhof(q, lambda, h = hr, delta = delta2, 
        epsabs = 1e-10, epsrel = 1e-8, limit = 1e6)$Qq)

ccdf_imhof4 <- sapply(q_vals, function(q) 
  imhof(q, lambda, h = hr, delta = delta2, 
        epsabs = 1e-11, epsrel = 1e-8, limit = 1e6)$Qq)


ccdf_farebrother <- sapply(q_vals, function(q) 
  farebrother(q, lambda, h = hr, delta = delta2, 
              maxit = 1e5, eps = 1e-12, mode = 1)$Qq)



###############################################################################
# Plotting code for the CCDF methods
###############################################################################
  
plot(q_vals, ccdf_imhof1, type = "o", 
     ylim = c(1e-8,1e-4), 
     #xlim = c(0,120), 
     xlab = "q", 
     ylab = "CCDF", 
     lwd = 2, 
     col = "black",
     log = "y",
     main = paste("Comparison of CCDF  - Examples of Imhof's Failure", exNumber))

lines(q_vals, ccdf_imhof2,       col = "red",   lwd = 2, lty = 2)
lines(q_vals, ccdf_imhof3,          col = "blue",  lwd = 1, lty = 3,pch = 20, type = "o")
lines(q_vals, ccdf_imhof4,          col = "purple",lwd = 2, lty = 5)
lines(q_vals, ccdf_farebrother,  col = "green", lwd = 3, lty = 4)

# Add grid for better visualization on log scale
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")

# # Only plot Zhang, LPB, WF if they are not NA
# if(!all(is.na(ccdf_zhang))) {
#   lines(q_vals, ccdf_zhang,  col = "orange", lwd = 2, lty = 6)
# }

legend("topleft", 
       legend = c("Imhof (epsabs = 1e-8)", "Imhof (epsabs = 1e-9)", 
                  "Imhof (epsabs = 1e-10)", "Imhof (epsabs = 1e-11)", "Farebrother"),
       col = c("black", "red", "blue", 
               "purple", "green" ),
       lty = c(1,2,3,4,5),
       bty = "n", cex = 0.8)


