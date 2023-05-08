
#' pfcre
#'
#' \code{pfcre} implements the Penalized Flexible Correlated Random Effects estimator
#'
#' @param formula an object of class \code{\link{formula}}: "y~x1+x2+x3..."
#' @param addcovars an object of class \code{\link{formula}}: "~z1+z2+z3..."
#' @param data an object of class \code{\link{data.frame}} containing the variables to be
#' used in the model, including group ids
#' @param id a variable in \code{data} that identifies groups
#' @param family a descriptions of the error distribution and link function to be used in the
#' model. See \code{\link{family}} for details of family functions.
#' @param lambda is the parameter value for the \code{scadpen} penalty. "a" is set to 3.7
#' @param weights an optional vector of prior weights to be used in the fitting process.
#' @param initial.beta starting values for the parameters in the linear predictor (optional)
#' @param mustart starting values for the vector of means (optional)
#' @param eta.new starting values for the linear predictor (optional)
#' @return Returns an object with the following elements:
#' \item{beta}{ a named vector of estimated main coefficients}
#' \item{gamma}{ a named vector of estimated coefficients for the polynomial}
#' \item{sdeta}{ the variance of the random effect}
#' \item{sdbeta}{ a named vector of the standard error of the \code{beta} coefficients}
#' \item{sdgamma}{ a named vector of the standard error of the \code{gamma} coefficients}
#' \item{vcov}{ a matrix of the covariance of \code{beta} and \code{gamma}}
#' \item{nobs}{ the number of observations}
#' \item{loglik}{ the unpenalized log-likelihood}
#' \item{formula}{ the formula used in the call}
pfcre <- function (formula = NULL, addcovars = NULL, data = NULL, id = NULL,
                   family = NULL, intercept = TRUE, weights = rep(1, nobs), initial.beta,
                   mustart, eta.new, lambda = NULL, max.steps = 1000, degree = 2){
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("lme4 needed for this function to work. Please install it.",
  call. = FALSE)
}
if (!requireNamespace("optimx", quietly = TRUE)) {
    stop("lme4 needed for this function to work. Please install it.",
         call. = FALSE)
  }

# Set matrices and elements
data.prep <- dataprep(formula, addcovars, data, "id", degree = degree)
nobs <- nrow(data)
y <- data.prep$y
x <- data.prep$xx
k <- data.prep$k
depvar <- data.prep$depvar

# Estimate model
PFCREfit <- pfcre1(y = y, x = x, k = k, id = id, family = family, lambda = lambda,
                 depvar = depvar, max.steps = max.steps)

# Extract elements from estimates
finform <- PFCREfit$finform
depvar <- PFCREfit$depvar
vars <- PFCREfit$vars
x <- PFCREfit$x
k <- PFCREfit$k

# Optimization options
nlopt <- function(par, fn, lower, upper, control) {
  .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper,
                            opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                        maxeval = 5000, xtol_abs = 1e-6, ftol_abs = 1e-6))
  list(par = res$solution,
       fval = res$objective,
       conv = if (res$status > 0) 0 else res$status,
       message = res$message
  )
}

# Final fit with selected variables
cat(".")
hybrid <- suppressWarnings(lme4::glmer(finform, data = data.frame(data[,c(depvar,id)], x[,vars]),
                                       family = family,
                                       control = lme4::glmerControl(calc.derivs = FALSE,optimizer = "nloptwrap")))

# Extract final estimation results
beta <- summary(hybrid)$coefficients[1:(k+1),1]
if (length(summary(hybrid)$coefficients[,1])<k+2){
  gamma <- NA
} else{
  gamma <- summary(hybrid)$coefficients[(k+2):dim(summary(hybrid)$coefficients)[1],1]
  gamma <- as.vector(gamma)
  names(gamma) <- rownames(summary(hybrid)$coefficients)[(k+2):dim(summary(hybrid)$coefficients)[1]]
}

sdbeta <- summary(hybrid)$coefficients[1:(k+1),2]
if (length(summary(hybrid)$coefficients[,1])<k+2){
  sdgamma <- NA
} else {
  sdgamma <- summary(hybrid)$coefficients[(k+2):dim(summary(hybrid)$coefficients)[1],2]
}

sdeta <- hybrid@theta

if (length(summary(hybrid)$coefficients[,1])<k+2){
  vcov <- matrix(vcov(summary(hybrid))@x,(length(beta)))
}
else {
  vcov <- matrix(vcov(summary(hybrid))@x, (length(beta)+length(gamma)))
}

loglik <- logLik(hybrid)
nobs <- dim(hybrid@frame)[1]
groups <- length(unique(hybrid@frame$id))

allform <- list(formula, addcovars)

fitted<-list(beta = beta, gamma = gamma, sdeta = sdeta, sdbeta = sdbeta, sdgamma = sdgamma,
             vcov = vcov, nobs = nobs, groups = groups, loglik = loglik, formula = allform )
fitted
}
