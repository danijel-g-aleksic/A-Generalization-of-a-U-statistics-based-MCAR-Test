# function that calculates Tn on one pair
TnXR <- function(X, R)
{
  XR <- cbind(X, R)
  n <- length(X)
  Tn_tilde <- mean(X)*mean(R) - mean(X*R)
  Tn <- n*Tn_tilde/(n-1)
  
  return(Tn)
}


# Function that returns An
An <- function(XY) {
  
  colind_x <- c()
  for (i in 1:ncol(XY)) {
    if(!any(is.na(XY[, i]))) {
      colind_x <- c(colind_x, i)
    }
  }
  
  ## FORMAT CHECK ###
  
  
  # wee need to have at least one incomplete column
  if(length(colind_x) == ncol(XY)) {
    stop("There is no missing data.")
  }
  
  X <- XY[, colind_x]
  R <- 1*(!is.na(XY[, -colind_x]))
  p <- length(colind_x)
  q <- ncol(XY) - p
  n <- nrow(XY)
  
  
  
  vec <- rep(0, p*q)
  mask <- cbind(vec, vec)
  
  # we need to cover special cases of p = 1 and q = 1
  if (p == 1) {
    for (k in 1:(p*q)) {
      vec[k] <- TnXR(X, R[, k])
    }
    S <- cov(R) * var(X)
    An <- vec %*% inv(S) %*% vec
    An <- n * An
    return(as.numeric(An))
  }
  
  if (q == 1) {
    for (k in 1:(p*q)) {
      vec[k] <- TnXR(X[, k], R)
    }
    S <- cov(X) * var(R)
    An <- vec %*% inv(S) %*% vec
    An <- n * An
    return(as.numeric(An))
  }
  
  k <- 1
  for(i in 1:p) {
    for (j in 1:q) {
      vec[k] <- TnXR(X[, i], R[, j])
      mask[k, ] <- c(i, j)
      k <- k + 1
    }
  }
  
  CVX <- cov(X)
  CVR <- cov(R)
  
  S <- matrix(rep(0, p*q*p*q), nrow = p*q)
  for(i in 1:(p*q)) {
    for(j in 1:(p*q)) {
      S[i, j] <- CVX[mask[i, 1], mask[j, 1]] * CVR[mask[i, 2], mask[j, 2]]
    }
  }
  
  An <- vec %*% pseudoinverse(S) %*% vec
  An <- as.numeric(n * An)
  pval <- pchisq(An, df = p*q, lower.tail = FALSE)
  return(pval)
}
