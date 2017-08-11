#'
#' @title Checks what values are not missing
#' @description This fucntion is a wrapper for: !is.na()
#' @details if the length of input vector is less than the set filter
#' a missing value is returned.
#' @param x atomic vector, list and pairlist to be tested
#' @return a suitable index vector for use with x.
#' @export
#'
notNaDS <- function (x) {
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  #check <- dsBase::isValidDS(x)
  
  # return missing value if the input vector is not valid
  if(TRUE){
    result <- !is.na(x)
  }else{
    result <- NA
  }
  
  return(result)
}
