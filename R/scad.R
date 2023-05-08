
scadpen <- function (lambda, beta){
  lambda[2] <- 3.7
  names(lambda) <- c("lambda", "a")
  first.derivative <- function(beta) {
    if (is.null(beta))
      stop("'beta' must be the current coefficient vector \n")
    lambda1 <- lambda[1]
    a <- lambda[2]
    theta <- abs(beta)
    p <- length(beta)
    help1 <- sapply(theta, function(theta) {
      max(a * lambda1 - theta, 0)/((a - 1) * lambda1)
    })
    lambda1 * ((theta <= lambda1) + help1 * (theta > lambda1))
  }
  getmat <- function(beta = beta, c1 = 1e-08) {
    if (is.null(beta))
      stop("'beta' must be the current coefficient vector \n")
    if (c1 < 0)
      stop("'c1' must be non-negative \n")
    penmat <- diag(first.derivative(beta)/(sqrt(beta^2 +c1) + 1e-06), length(beta))
  }
  getmat(beta)
}
