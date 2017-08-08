#'
#' @title Negate
#' @description Negate
#' @details if the length of input vector is less than the set filter
#' a missing value is returned.
#' @param x
#' @return negation
#' @export
#'
notNA <- function (x) {
  
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  #check <- isValidDS(xvect)
  
  # return missing value if the input vector is not valid
  #if(check){
    result <- !is.na(x)
  #}else{
  #  result <- NA
  #}
  
  return(result)
}
