
pfcre1 <- function(y = NULL, x = NULL, k = NULL, depvar = NULL, id = NULL ,
                   family = NULL, intercept = TRUE, weights = rep(1, nobs), initial.beta,
                   mustart, eta.new, lambda=NULL, max.steps = NULL)
{


# Convergence Criteria and Initialization
  gamma <- 1
  x <- as.matrix(x)
  converged <- FALSE
  n.iter <- max.steps
  eps <- 1.5e-06
  stop.at <- n.iter
  p <- ncol(x)
  nobs <- nrow(x)
  converged <- FALSE
  glm1 <- glm.fit(x,y,family=binomial(link = "logit"))
  beta.mat <- matrix(0, nrow = n.iter, ncol = p)
  if (missing(initial.beta))
    initial.beta <- glm1$coefficients
  else eta.new <- drop(x %*% initial.beta)
  if (missing(mustart)) {
    etastart <- drop(x %*% initial.beta)
    eval(family$initialize)
  }
  if (missing(eta.new))
    eta.new <- family$linkfun(mustart)

# Algorithm iterations
    for (i in 1:n.iter) {
    beta.mat[i, ] <- initial.beta
    mu.new <- family$linkinv(eta.new)
    d.new <- family$mu.eta(eta.new)
    v.new <- family$variance(mu.new)
    weights <- d.new/sqrt(v.new)
    x.star <- weights * x
    y.tilde.star <- weights * (eta.new + (y - mu.new)/d.new)

    #######################################################
    # Penalty matrix
    A.lambda <- scadpen(lambda = lambda, beta = initial.beta)
    A.lambda[1:(k+1),1:(k+1)]=0
    p.imat.new <- crossprod(x.star) + A.lambda
    chol.pimat.new <- chol(p.imat.new)
    inv.pimat.new <- chol2inv(chol.pimat.new)
    beta.new <- gamma * drop(inv.pimat.new %*% t(x.star) %*%
                               y.tilde.star) + (1 - gamma) * beta.mat[i, ]
    if ((sum(abs(beta.new - initial.beta))/sum(abs(initial.beta)) <=
         eps)) {
      converged <- TRUE
      stop.at <- i
      if (i < n.iter)
        break
    }
    else {
      initial.beta <- beta.new
      eta.new <- drop(x %*% beta.new)
    }
    }

# Res. and fit eval.
  Hatmat <- x.star %*% inv.pimat.new %*% t(x.star)
  tr.H <- sum(diag(Hatmat))
  dev.m <- sum(family$dev.resids(y, mu.new, weights))
  aic.vec <- dev.m + 2 * tr.H
  bic.vec <- dev.m + log(nobs) * tr.H

# Variable selection & Outputs
  coefficients <- beta.new
  vars <- which(abs(beta.new)>5e-03)[-1]

  finform <- as.formula(paste(depvar,"~", paste(colnames(x[,vars]), collapse="+"), "+(1|","id",")"))

  fit <- list(coefficients = beta.new, beta.mat = beta.mat[1:stop.at, ], tr.H = tr.H, fitted.values = mu.new,
              family = family, Amat = A.lambda, converged = converged, stop.at = stop.at,
              m.stop = stop.at, linear.predictors = eta.new, weights = weights^2,
              p.imat = p.imat.new, inv.pimat = inv.pimat.new, x.star = x.star,
              v.new = v.new)
  PFCREfit <- list(finform = finform, depvar = depvar, x = x, vars = vars, id = id,
                   data = data, family = family, k = k, fit = fit)
  PFCREfit
}
