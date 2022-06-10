# Functions for model simulation, analysis, and visualization
# zachmiller@uchicago.edu

library("deSolve")


sample_P <- function(n, check.feas = TRUE, check.stab = TRUE, symmetric = FALSE, 
                     m = 1, distribution = runif, diag.distribution = distribution, ...) {
  
  # Generate a P matrix at random, and possibly according to some constraints
  # 
  # Input
  # n : number of species
  # check.feas : should the output satisfy P^{-1} m > 0 and m 1^T P^{-1} 1 < 1 ? (keep sampling until found)
  # check.stab : should the output have a stable coexistence equilibrium? (keep sampling until found)
  #  (NOTE: this option is currently only available for symmetric matrices)
  # symmetric : should the output be symmetric?
  # m : one of (i) a vector of length n or (ii) a single value to be used for all species
  # distribution : a function to generate random deviates, from which elements of P are sampled i.i.d.
  #     (subject to symmetry or separate distribution for diagonal)
  # diag.distribution : a function to generate i.i.d. random deviates for the diagonal elements of P
  # ... : any additional parameters for distribution functions
  #
  # Output
  # an n x n matrix P satisfying input constraints, and logicals feas (does the output have feasible equilibrium)
  #     and stab (does the output have stable equilibrium)
  
  success <- FALSE
  while(!success) { # sample until matrix with desired properties is found
    
    stab <- FALSE
    feas <- FALSE
    
    P <- matrix(distribution(n = n*n, ...), n, n) # draw entries i.i.d. from distribution function
    if (symmetric) P <- (P + t(P)) / 2 # enforce symmetry if required
    diag(P) <- diag.distribution(n, ...) # sample diagonal from another distribution (same as off-diagonal by default)
    
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


slow_system <- function(t, state, pars){
  
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
