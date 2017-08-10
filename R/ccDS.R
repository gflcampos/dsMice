#'
#'@title Select complete cases
#'@description Extracts the complete cases, also known as \emph{listwise deletion}.
#'\code{cc(x)} is similar to 
#'\code{na.omit(x)}, but returns an object of the same class 
#'as the input data. Dimensions are not dropped. For extracting
#'incomplete cases, use \code{\link{ici}}.
#'
#'@param x An \code{R} object. Methods are available for classes
#'\code{mids}, \code{data.frame} and \code{matrix}. Also, \code{x} 
#'could be a vector.
#'@return A \code{vector}, \code{matrix} or \code{data.frame} containing the data of the complete cases.
#'@author Stef van Buuren, 2017.
#'@seealso \code{\link{na.omit}}, \code{\link{cci}}, \code{\link{ici}}
#'@keywords univar
#'@examples
#'
#'# cc(nhanes)   # get the 13 complete cases
#'# cc(nhanes$bmi) # extract complete bmi
#'@export
ccDS <- function(x) {
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  #check <- isValidDS(x)
  
  if(TRUE) {
    result <- mice::cc(x)
  } else { # return NA if the input vector is not valid
    result <- NA
  }
  
  return(result)
}
