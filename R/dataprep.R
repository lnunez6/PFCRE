
dataprep <- function(formula = NULL, addcovars = NULL, data = NULL, id = NULL, degree = NULL){

  y <- model.frame(formula, data)[,1]
  x <- model.frame(formula, data)[,-1]
  depvar <- colnames(model.frame(formula, data))[1]
  xb <- apply(x, 2 , function(i){ave(i, data[ , id])})
  colnames(xb) <- paste("A.",colnames(xb),sep = "")

  if(is.null(addcovars)){ xbz <- xb
  }
  else {xbz <- cbind(xb,model.frame(addcovars,data))
  }
  xbz <- data.frame(xbz)


  Z1 <- model.matrix(as.formula(paste("~(", paste(colnames(xbz), collapse = "+"),")^2")),xbz)[,-1]
  Z2 <- outer(as.matrix(xbz), 2, "^")[,, 1]
  colnames(Z2) <- paste("SQ", colnames(Z2), sep = "")

  Z <- cbind(Z1, Z2)

  # Remove collinear
  qrZ <-  qr(Z, tol = 1e-8, LAPACK = FALSE)
  rankZ <- qrZ$rank
  keep <- qrZ$pivot[seq_len(rankZ)]
  Z <- Z[, keep]

  colnames(Z) <- paste("pol", 1:dim(Z)[2], sep = "")

  xx <- cbind(1, x, Z)
  colnames(xx)[1] <- "intercept"
  k <- dim(x)[2]

  data.prep <- list(y = y, xx = xx, depvar = depvar, k = k)
}


