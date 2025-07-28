#############################
## Helper Functions
#############################

# CGF: K_Q(t)
KQ <- function(t, lambda, nu, h2) {
  # K_Q(t) = -1/2 * sum_{i} nu_i * ln(1 - 2*lambda_i*t)
  #         + t * sum_{i} [ h2_i * lambda_i / (1 - 2*lambda_i*t ) ]
  part1 <- -0.5 * sum( nu * log(1 - 2*lambda*t) )
  part2 <- t * sum( h2 * lambda / (1 - 2*lambda*t) )
  part1 + part2
}

# 1st derivative of K_Q: K_Q'(t)
KQ_prime <- function(t, lambda, nu, h2) {
  # From the derived formula:
  #   K_Q'(t) = sum_{i} [ nu_i * lambda_i / (1 - 2*lambda_i*t ) ]
  #           + sum_{i} [ h2_i * lambda_i / (1 - 2*lambda_i*t ) ]
  #           + 2 * t * sum_{i} [ h2_i * lambda_i^2 / (1 - 2*lambda_i*t)^2 ]
  # Combine the first two sums:
  sum1 <- sum( (nu + h2) * lambda / (1 - 2*lambda*t) )
  sum2 <- 2 * t * sum( h2 * lambda^2 / (1 - 2*lambda*t)^2 )
  sum1 + sum2
}

# 2nd derivative of K_Q: K_Q''(t)
KQ_doubleprime <- function(t, lambda, nu, h2) {
  # From the derived formula:
  #   K_Q''(t) = 2 * sum_{i} [ lambda_i^2 * (nu_i + 2*h2_i) / (1 - 2*lambda_i*t)^2 ]
  #            + 8 * t * sum_{i} [ h2_i * lambda_i^3 / (1 - 2*lambda_i*t)^3 ]
  partA <- 2 * sum( lambda^2 * (nu + 2*h2) / (1 - 2*lambda*t)^2 )
  partB <- 8 * t * sum( h2 * lambda^3 / (1 - 2*lambda*t)^3 )
  partA + partB
}

# 3rd derivative of K_Q: K_Q'''(t)
KQ_tripleprime <- function(t, lambda, nu, h2) {
  # From the derived formula:
  #  K_Q'''(t) = 8 * sum_{i} [ lambda_i^3 (nu_i + 3*h2_i) / (1 - 2*lambda_i*t)^3 ]
  #              + 48 * t * sum_{i} [ h2_i * lambda_i^4 / (1 - 2*lambda_i*t)^4 ]
  partA <- 8 * sum( lambda^3 * (nu + 3*h2) / (1 - 2*lambda*t)^3 )
  partB <- 48 * t * sum( h2 * lambda^4 / (1 - 2*lambda*t)^4 )
  partA + partB
}

#############################
## Main Function: sap()
#############################

spa <- function(q,lambda, nu, h2,eps=1e-3) {
  # 1) Basic checks
  stopifnot(length(lambda) == length(nu), length(nu) == length(h2))
  N <- length(lambda)
  if (N < 1) stop("Vectors 'lambda', 'nu', 'h2' must be non-empty.")
  
  # 2) Compute the mean of Q: mu = sum_i [ lambda_i * (nu_i + h2_i) ]
  q_mean <- sum(lambda * (nu + h2))
  
  ########################################
  ## Case A: q is essentially the mean
  ########################################
  if (abs(q - q_mean) < eps) {
    # Evaluate  K''(0) and K'''(0):
    K2_0 <- KQ_doubleprime(0, lambda, nu, h2)
    K3_0 <- KQ_tripleprime(0, lambda, nu, h2)
    
    # F_Q(mu) ~ 1/2 + [K_Q^{(3)}(0)] / [6 * sqrt(2*pi) * (K_Q^{(2)}(0))^(3/2)]
    val <- 0.5 + (K3_0) / (6 * sqrt(2*pi) * (K2_0)^(3/2))
    return(val)
    
  } else {
    ########################################
    ## Case B: solve K'(t) = q  for t
    ########################################
    
    # 1) Build a function for uniroot:
    eqfun <- function(z) {
      KQ_prime(z, lambda, nu, h2) - q
    }
    
    # 2) Decide the search interval based on sign of lambda:
    lmax <- max(lambda)
    lmin <- min(lambda)
    
    if (all(lambda > 0)) {
      # Case 1: All lambdas positive
      lowerB <- -1e10
      upperB <- 0.99999 / (2 * lmax)
    } else if (all(lambda < 0)) {
      # Case 2: All lambdas negative
      lowerB <- 0.99999 / (2 * lmin)
      upperB <- 1e10
    } else {
      # Case 3: Mixed signs
      lowerB <- 0.99999 / (2 * lmin)
      upperB <- 0.99999 / (2 * lmax)
    }
    
    # 3) Solve for zeta1
    zeta1 <- uniroot(eqfun, lower=lowerB, upper=upperB, tol=1e-8)$root
    
    ########################################
    ## Evaluate the saddlepoint approximation
    ########################################
    
    # w(zeta) = sign(zeta)* sqrt( 2*(zeta*q - KQ(zeta)) )
    w_zeta <- function(z) {
      val <- z*q - KQ(z, lambda, nu, h2)
      # sign(z)*sqrt(...) is the usual formula
      sign(z) * sqrt(2 * val)
    }
    
    # v(zeta) = z * sqrt( KQ''(z) )
    v_zeta <- function(z) {
      z * sqrt(KQ_doubleprime(z, lambda, nu, h2))
    }
    
    wVal <- w_zeta(zeta1)
    vVal <- v_zeta(zeta1)
    
    # The final CDF ~ Phi( w + (1/w)* log(v/w) )
    argument <- wVal + (1 / wVal) * log(vVal / wVal)
    cdf_approx <- pnorm(argument, mean=0, sd=1)
    return(cdf_approx)
  }
  
}

