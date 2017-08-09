#'
#'@title placeholder
#'@description placeholder
#'@param y Vector to be imputed
#'@param ry Logical vector of length \code{length(y)} indicating the 
#'the subset \code{y[ry]} of elements in \code{y} to which the imputation 
#'model is fitted. The \code{ry} generally distinguishes the observed 
#'(\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
#'@param x Numeric design matrix with \code{length(y)} rows with predictors for 
#'\code{y}. Matrix \code{x} may have no missing values.
#'@param wy Logical vector of length \code{length(y)}. A \code{TRUE} value 
#'indicates locations in \code{y} for which imputations are created.
#'@param donors The size of the donor pool among which a draw is made. 
#'The default is \code{donors = 5L}. Setting \code{donors = 1L} always selects 
#'the closest match, but is not recommended. Values between 3L and 10L 
#'provide the best results in most cases (Morris et al, 2015).
#'@param matchtype Type of matching distance. The default choice 
#'(\code{matchtype = 1L}) calculates the distance between 
#'the \emph{predicted} value of \code{yobs} and 
#'the \emph{drawn} values of \code{ymis} (called type-1 matching). 
#'Other choices are \code{matchtype = 0L} 
#'(distance between predicted values) and \code{matchtype = 2L} 
#'(distance between drawn values).
#'@param ridge The ridge penalty used in \code{.norm.draw()} to prevent 
#'problems with multicollinearity. The default is \code{ridge = 1e-05}, 
#'which means that 0.01 percent of the diagonal is added to the cross-product. 
#'Larger ridges may result in more biased estimates. For highly noisy data 
#'(e.g. many junk variables), set \code{ridge = 1e-06} or even lower to 
#'reduce bias. For highly collinear data, set \code{ridge = 1e-04} or higher.
#'@param \dots Other named arguments.
#'@return Vector with imputed data, same type as \code{y}, and of length 
#'\code{sum(wy)}
#'@details details
#'@author Stef van Buuren, Karin Groothuis-Oudshoorn
#'@return a numeric, the statistical mean
#'@export
#'
pmmDS <- function (y, ry, x, wy = NULL, donors = 5, 
                   matchtype = 1L, ridge = 1e-05, ...)
{
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  #check <- isValidDS(xvect)
  check <- TRUE
  
  # return NA if the input vector is not valid
  if(check){ # disable datashield check for now
    if (is.null(wy))
      wy <- !ry
    x <- cbind(1, as.matrix(x))
    ynum <- y
    if (is.factor(y)) 
      ynum <- as.integer(y)
    parm <- norm.drawDS(ynum, ry, x, ridge = ridge, ...)
    
    if (matchtype == 0L) {
      yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
      yhatmis <- x[wy, , drop = FALSE] %*% parm$coef
    }
    if (matchtype == 1L) {
      yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
      yhatmis <- x[wy, , drop = FALSE] %*% parm$beta
    }
    if (matchtype == 2L) {
      yhatobs <- x[ry, , drop = FALSE] %*% parm$beta
      yhatmis <- x[wy, , drop = FALSE] %*% parm$beta
    }
    #idx <- matcher(yhatobs, yhatmis, k = donors)
    
    result <- yhatmis#y[ry][idx]
  }else{
    result <- NA
  }
  
  return(result)
}
