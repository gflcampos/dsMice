#' Draws values of beta and sigma by Bayesian linear regression
#' 
#' This function draws random values of beta and sigma under the Bayesian 
#' linear regression model as described in Rubin (1987, p. 167). This function
#' can be called by user-specified imputation functions.
#' 
#'@param y Incomplete data vector of length \code{n}
#'@param ry Vector of missing data pattern (\code{FALSE}=missing,
#'\code{TRUE}=observed)
#'@param x Matrix (\code{n} x \code{p}) of complete covariates.
#'@param ridge A small numerical value specifying the size of the ridge used. 
#' The default value \code{ridge = 1e-05} represents a compromise between stability
#' and unbiasedness. Decrease \code{ridge} if the data contain many junk variables.
#' Increase \code{ridge} for highly collinear data. 
#'@param ... Other named arguments.
#'@return A \code{list} containing components \code{coef} (least squares estimate),
#'\code{beta} (drawn regression weights) and \code{sigma} (drawn value of the 
#'residual standard deviation).
#'@references
#'Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York: Wiley.
#'@author Stef van Buuren, Karin Groothuis-Oudshoorn, 2000
#'@export
norm.drawDS <- function(y, ry, x, ridge = 1e-05, ...) {
  xobs <- x[ry, ]
  yobs <- y[ry]
  xtx <- crossprod(xobs)
  pen <- ridge * diag(xtx)
  if (length(pen) == 1)
    pen <- matrix(pen)
  v <- solve(xtx + diag(pen))
  coef <- t(yobs %*% xobs %*% v)
  residuals <- yobs - xobs %*% coef
  df <- max(sum(ry) - ncol(x), 1)
  sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))
  #beta.star <- coef + (t(chol(sym(v))) %*% rnorm(ncol(x))) * sigma.star #FAILS!
  #parm <- list(coef, beta.star, sigma.star)
  #names(parm) <- c("coef", "beta", "sigma")
  return(coef)
}