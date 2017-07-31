#'
#' @title Computes statistical mean of a vector
#' @description Calculates the mean value.
#' @details if the length of input vector is less than the set filter
#' a missing value is returned.
#' @param xvect a vector
#' @return a numeric, the statistical mean
#' @author Gaye, A.
#' @export
#'
pmmDS <- function (y, ry, x, wy = NULL, donors = 5, 
                   matchtype = 1L, ridge = 1e-05, ...)
{
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  #check <- isValidDS(xvect)
  check <- TRUE
  
  # return NA if the input vector is not valid
  if(check){
    if (is.null(wy)) 
      wy <- !ry
    x <- cbind(1, as.matrix(x))
    ynum <- y
    if (is.factor(y)) 
      ynum <- as.integer(y)
    parm <- norm.drawDS(ynum, ry, x, ridge = ridge, ...)
    
    result <- parm
  }else{
    result <- NA
  }
  
  return(result)
}
