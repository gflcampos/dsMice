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
    
    sourceCpp(code='
        #include <Rcpp.h>
        #include <algorithm>
        using namespace std;
        using namespace Rcpp;
        
        // [[Rcpp::export]]
        IntegerVector matcher(NumericVector obs, NumericVector mis, int k) {
          // fast predictive mean matching algorithm
          // for each of the n0 elements in mis
          // 1) calculate the difference with obs
          // 2) add small noise to break ties
          // 3) find the k indices of the k closest predictors
          // 4) randomly draw one index
          // and return the vector of n0 matched positions 
          // SvB 26/01/2014
          
          // declarations
          int jj;
          int n1 = obs.size();
          int n0 = mis.size();
          double dk = 0; 
          int count = 0;
          int goal = 0;
          NumericVector d(n1);
          NumericVector d2(n1);
          IntegerVector matched(n0);
          
          // restrict 1 <= k <= n1
          k = (k <= n1) ? k : n1;
          k = (k >= 1) ? k : 1;
          
          // in advance, uniform sample from k potential donors
          NumericVector which = floor(runif(n0, 1, k + 1));
          NumericVector mm = range(obs);
          double small = (mm[1] - mm[0]) / 65536;
          
          // loop over the missing values
          for(int i = 0; i < n0; i++) {
          
          // calculate the distance and add noise to break ties
          d2 = runif(n1, 0, small);
          dk = mis[i];
          for (int j = 0; j < n1; j++) d[j] = std::abs(obs[j] - dk) + d2[j];
          
            // find the kth lowest value in d
            for (int j = 0; j < n1; j++) d2[j] = d[j];
            std::nth_element (d2.begin(), d2.begin() + k - 1, d2.end());
            
            // find index of donor which[i]
            dk = d2[k-1];
            count = 0;
            goal = (int) which[i];
            for (jj = 0; jj < n1; jj++) {
              if (d[jj] <= dk) count++;
              if (count == goal) break;
            }
            
            // and store the result
            matched[i] = jj;
          }
    
          // increase index to offset 1
          return matched + 1;
        }
        
        static R_CallMethodDef callMethods[]  = {
          {"matcher", (DL_FUNC) &matcher, 3},
          {NULL, NULL, 0}
        };
        
        void attribute_visible R_init_mice(DllInfo *dll)
        {
          R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
          R_useDynamicSymbols(dll, FALSE);
          R_forceSymbols(dll, TRUE);
        }
    ')
    
    idx <- matcher(yhatobs, yhatmis, k = 5)
    
    result <- y[ry][idx]
  }else{
    result <- NA
  }
  
  return(result)
}
