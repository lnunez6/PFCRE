
#' pfcrems
#'
#' \code{pfcrems} performs model penalty selection for the Penalized Flexible Correlated Random Effects estimator
#'
#' @param formula an object of class \code{\link{formula}}: "y~x1+x2+x3..."
#' @param addcovars an object of class \code{\link{formula}}: "~z1+z2+z3..."
#' @param data an object of class \code{\link{data.frame}} containing the variables to be
#' used in the model, including group ids
#' @param id a variable in \code{data} that identifies groups
#' @param family a descriptions of the error distribution and link function to be used in the
#' model. See \code{\link{family}} for details of family functions.
#' @param lambda.cand is a vector of candidate values for the \code{scadpen} penalty. "a" is set to 3.7
#' @param weights vector of prior weights to be used in the fitting process (optional)
#' @param initial.beta starting values for the parameters in the linear predictor (optional)
#' @param mustart starting values for the vector of means (optional)
#' @param eta.new starting values for the linear predictor (optional)
#' @return Returns an object with the following elements:
#' \item{Model}{ a list with the estimates from the best model selected via Akaike Information Criterion}
#' \item{lambda.cand}{ a vector with the original candidate values for the \code{scadpen} penalty}
#' \item{cv}{ a vector with the model fit for the \code{labmda.cand} elements}
#' @export



pfcrems <- function(formula = NULL, addcovars = NULL, data = NULL, id = NULL, family = NULL,
                    intercept = TRUE, weights = rep(1, nobs), initial.beta,
                    mustart, eta.new,lambda.cand = NULL, max.steps = 1000, degree = 2)
  {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("lme needed for this function to work. Please install it.",
           call. = FALSE)
    }
# Data preliminaries
    data.prep <- dataprep(formula, addcovars, data, id, degree = degree)
    nobs <- nrow(data)

    #############################################################################################
    #############################################################################################
    fit_model <- list()
    i = 1
    for (i in 1:length(lambda.cand)){
      fit_model[[i]] <- pfcre(formula = formula, addcovars = addcovars,
                              data = data, id = id, family=family,
                              intercept = intercept, weights = weights,
                              initial.beta = init.beta, mustart, eta.new, lambda = lambda.cand[i],
                              max.steps=max.steps,
                              degree=degree)
      init.beta <- c(fit_model[[i]]$beta,fit_model[[i]]$gamma)
      i = i+1
    }

    fit_aic <- unlist(lapply(fit_model, function(i) -2*i$loglik+2*(length(i$gamma)+length(i$beta)+1)))
    fit_bic <- unlist(lapply(fit_model, function(i) -2*i$loglik+log(nobs)*(length(i$gamma)+length(i$beta)+1)))
    selected <- which(fit_aic == min(fit_aic))[1]

    pf <- fit_model[[selected]]
    #############################################################################################
    #############################################################################################

    lo <- lambda.cand[selected]
    if (lo == min(lambda.cand)) warning("All lambda candidates too high, choose smaller values \n",
                                        call. = FALSE)
    if (lo == max(lambda.cand)) warning("All lambda candidates too low, choose higher values \n",
                                        call. = FALSE)

    # Best model fit
#    pf <- pfcre(formula = formula, addcovars = addcovars, data = data, id = id, family = family,
#                intercept = intercept, weights = weights, initial.beta, mustart, eta.new, lambda = lo,
#                max.steps = max.steps, degree = degree)

    # Final output
    pfcrems <- list(Model = pf, lambda.cand = lambda.cand, cv = fit_aic)
    pfcrems
  }
