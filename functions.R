# Functions for model simulation, analysis, and visualization
# zachmiller@uchicago.edu

library("deSolve")


sample_P <- function(n, check.feas = TRUE, check.stab = TRUE, symmetric = FALSE, 
                     m = 1, eps = 0.5, distribution = runif, diag.distribution = distribution, ...) {
  
  # Flexible function to generate a P matrix at random, and possibly according to some constraints
  # 
  # Input
  # n : number of species
  # check.feas : should the output satisfy P^{-1} m > 0 and m 1^T P^{-1} 1 < 1 ? (keep sampling until found)
  # check.stab : should the output have a stable coexistence equilibrium? (keep sampling until found)
  #  (NOTE: this option is currently only available for symmetric matrices)
  # symmetric : should the output be symmetric?
  # m : one of (i) a vector of length n or (ii) a single value to be used for all species or (iii) the string "proportion" (see eps)
  # eps : if m is chosen as a proportion of m_max, eps sets the proportion (e.g. eps = 0.1 means m = 10% of m_max)
  # distribution : a function to generate random deviates, from which elements of P are sampled i.i.d.
  #     (subject to symmetry or separate distribution for diagonal)
  # diag.distribution : a function to generate i.i.d. random deviates for the diagonal elements of P
  # ... : any additional parameters for distribution functions
  #
  # Output
  # an n x n matrix P satisfying input constraints, and logicals feas (does the output have feasible equilibrium)
  #     and stab (does the output have stable equilibrium)
  
  use_proportions <- identical(m, "proportion") # check if using proportional m
  
  success <- FALSE
  while(!success) { # sample until matrix with desired properties is found
    
    stab <- FALSE
    feas <- FALSE
    
    P <- matrix(distribution(n = n*n, ...), n, n) # draw entries i.i.d. from distribution function
    if (symmetric) P <- (P + t(P)) / 2 # enforce symmetry if required
    diag(P) <- diag.distribution(n, ...) # sample diagonal from another distribution (same as off-diagonal by default)
    
    # if m = "proportion", compute max_m and then m as eps * max_m
    if (use_proportions) {
      max_m <- 1 / sum(rowSums(solve(P)))
      m <- eps * max_m
    }
    
    # if a single value of m is specified, repeat to obtain a vector of length n
    if (length(m) == 1) m <- rep(m, n)
    
    # check if P and m yield feasible equilibrium
    y <- solve(P, m)
    if (all(y > 0) & sum(y) < 1) feas <- TRUE
    
    # check if P and m yield stable equilibrium
    if (symmetric) { # for symmetric P, just check eigenvalue condition
      eP <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
      if (sum(Re(eP) > 0) == 1) stab <- TRUE # stable iff exactly one positive eigenvalue
    } else {
      if (check.stab) {
        warning("Stability is only checked for symmetric matrices. Returning P without checking stability.")
        stab <- TRUE
      }
    }
    
    # check that required properties are satisfied; if not, discard and try again
    if ((!check.feas | feas) & (!check.stab | stab)) success <- TRUE
  }
  
  return(list(P = P, feas = feas, stab = stab))
}


general_feedback_model <- function(t, state, pars){
  
  # dynamics given by Eq. 5 (for numerical integration)
  
  n <- pars$n
  P <- pars$P
  m <- pars$m
  d <- pars$d
  
  x <- state[1:n]
  y <- state[(n+1):(2*n)]
  w <- state[(2*n+1):(3*n)]
  
  dxdt <- x * (-m + P %*% y)
  dydt <- m * (w - y) - y * (t(P) %*% x)
  dwdt <- d * (x + y - w)
  
  return(list(c(dxdt, dydt, dwdt)))
}

unreduced_feedback_model <- function(t, state, pars){
  
  # dynamics given by Eq. S16 (for numerical integration)
  
  n <- pars$n
  P <- pars$P
  m <- pars$m
  d <- pars$d
  
  X <- matrix(state[1:n^2], n, n)
  y <- state[(n^2 + 1):(n^2 + n)]
  
  dXdt <- -m * X + diag(rowSums(X)) %*% P %*% diag(y) - d * X + diag(rowSums(d * X))
  dydt <- colSums(m * X) - diag(y) %*% t(P) %*% rowSums(X)
  
  return(list(c(as.vector(dXdt), dydt)))
}


slow_system <- function(t, state, pars){
  
  # dynamics given by Eq. 7 (for numerical integration)
  
  n <- pars$n
  P <- pars$P
  m <- pars$m
  d <- pars$d
  A <- pars$A
  b <- pars$b
  
  w <- state / sum(state)
  
  dwdt <- d * (A %*% w + b)
  
  return(list(c(dwdt)))
}
