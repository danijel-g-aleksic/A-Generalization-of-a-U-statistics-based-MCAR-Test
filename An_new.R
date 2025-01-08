# function that calculates Tn on one pair
TnXR_new <- function(X, R)
{
  XR <- cbind(X, R)
  XRo <- na.omit(XR)
  n <- nrow(XRo)
  Tn_tilde <- mean(XRo[, 1])*mean(XRo[, 2]) - mean(XRo[, 1] * XRo[, 2])
  Tn <- n*Tn_tilde/(n-1)
  
  return(Tn)
}

TnYR_new <- function(Y, R)
{
  n <- length(R)
  Yrepl <- Y
  Yrepl[is.na(Yrepl)] <- 0
  Tn_tilde <- mean(Yrepl)*mean(R) - mean(Yrepl*R)
  Tn <- n*Tn_tilde/(n-1)
  return(Tn)
}




# Function that returns An
An_new <- function(XY) {
  
  # indices of complete columns
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
  
  X <- data.frame(XY[, colind_x])
  Y <- XY[, -colind_x]
  R <- 1*(!is.na(XY[, -colind_x]))
  p <- length(colind_x)
  q <- ncol(XY) - p
  n <- nrow(XY)
  
  
  
  vec <- rep(0, p*q + q*(q-1))
  mask <- cbind(vec, vec)
  
  
  k <- 1
  for(i in 1:p) {
    for (j in 1:q) {
      vec[k] <- TnXR_new(X[, i], R[, j])
      mask[k, ] <- c(i, j)
      k <- k + 1
    }
  }
  
  for(i in 1:q) {
    for (j in 1:q) {
      if(i==j) {next}
      vec[k] <- TnYR_new(Y[, i], R[, j])
      mask[k, ] <- c(i, j)
      k <- k + 1
    }
  }
  
  CVX <- cov(X)
  CVR <- cov(R)

  
  SXR <- matrix(rep(0, p*q*p*q), nrow = p*q)
  for(i in 1:(p*q)) {
    for(j in 1:(p*q)) {
      SXR[i, j] <- CVX[mask[i, 1], mask[j, 1]] * CVR[mask[i, 2], mask[j, 2]]
    }
  }
  SYR <- matrix(rep(0, q*(q-1)*q*(q-1)), nrow = q*(q-1))
  for(i in 1:(q*(q-1))) {
    for(j in 1:(q*(q-1))) {
      Yu <- Y[, mask[i+p*q, 1]]
      Yr <- Y[, mask[j+p*q, 1]]
      Yu[is.na(Yu)] <- 0
      Yr[is.na(Yr)] <- 0
      SYR[i, j] <- cov(Yu, Yr) * CVR[mask[i+p*q, 2], mask[j+p*q, 2]] 
    }
  }
  
  S <- matrix(rep(0, (p*q + q*(q-1))*(p*q + q*(q-1))), nrow = p*q + q*(q-1))
  S[1:(p*q), 1:(p*q)] <- SXR
  S[(p*q + 1):(p*q + q*(q-1)), (p*q + 1):(p*q + q*(q-1))] <- SYR
  
  
  for(i in 1:(p*q)) {
    for(j in (p*q + 1):(p*q + q*(q-1))) {
      tmp <- cbind(X[, mask[i, 1]], Y[, mask[j, 1]])
      tmp <- na.omit(tmp)
      S[i, j] <- cov(tmp[, 1], tmp[, 2]) * CVR[mask[i, 2], mask[j, 2]] * mean(!is.na(Y[, mask[j, 1]])) # * CVR[mask[j, 2], mask[j, 2]]   
      S[j, i] <- S[i, j]
    }
  }
  
  An <- vec %*% pseudoinverse(S) %*% vec
  An <- as.numeric(n * An)
  pval <- pchisq(An, df = p*q + q*(q-1), lower.tail = FALSE)
  return(pval)
}
