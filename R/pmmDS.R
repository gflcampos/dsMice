#'
#'@title Calculates imputations for univariate missing data by predictive mean matching
#'@description This function performs imputation by predictive mean matching by executing the pmmDS function on the server-side.
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
#'@examples
#'
#'# In this example, we assume that the Opal server to which we are connecting, 
#'# has a table that contains the 'boys' data from the original mice package.
#'
#'# Load DataSHIELD libraries
#'library(dsBaseClient)
#'library(dsMiceClient)
#'
#'# Build login information
#'server <- c("server_name")
#'url <- c("opal_url")
#'user <- "username"
#'password <- "password"
#'table <- c("project_name.table_name")
#'logindata <- data.frame(server,url,user,password,table)
#'
#'# Login and assign the 'boys' dataset to varable 'D' on the server-side
#'opals <- datashield.login(logins=logindata, assign=TRUE)
#'
#'datashield.assign(opals, symbol="xname", value=as.symbol("c('age', 'hgt', 'wgt')"))
#'datashield.assign(opals, symbol="r", value=as.symbol("complete.cases(D[, xname])"))
#'datashield.assign(opals, symbol="x", value=as.symbol("D[r, xname]"))
#'datashield.assign(opals, symbol="y", value=as.symbol("D[r, 'tv']"))
#'datashield.assign(opals, symbol="ry", value=as.symbol("notNA(y)"))
#'
#'# Impute missing tv data
#'datashield.aggregate(opals, "pmmDS(y, ry, x)")
#'yimp <- ds.mice.pmm('y','ry','x')
#'length(yimp)
#'hist(yimp, xlab = 'Imputed missing tv')
#'
#'@export

pmmDS <- function (y, ry, x, wy = NULL, donors = 5, 
                   matchtype = 1L, ridge = 1e-05, ...) {
  # check if the input vector is valid (i.e. meets DataSHIELD privacy criteria)
  check <- dsBase:::isValidDS(y)
  
  if(check) {
    result <- mice::mice.impute.pmm(y, ry, x, wy=wy, donors=donors, matchtype=matchtype, ridge=ridge)
  } else { # return NA if the input is not valid
    result <- NA
  }
  
  return(result)
}
