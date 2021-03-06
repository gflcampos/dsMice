% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccDS.R
\name{ccDS}
\alias{ccDS}
\title{Select complete cases}
\usage{
ccDS(x)
}
\arguments{
\item{x}{An \code{R} object. Methods are available for classes
\code{mids}, \code{data.frame} and \code{matrix}. Also, \code{x} 
could be a vector.}
}
\value{
A \code{vector}, \code{matrix} or \code{data.frame} containing the data of the complete cases.
}
\description{
Extracts the complete cases, also known as \emph{listwise deletion}.
\code{cc(x)} is similar to 
\code{na.omit(x)}, but returns an object of the same class 
as the input data. Dimensions are not dropped. For extracting
incomplete cases, use \code{\link{ici}}.
}
\examples{

# cc(nhanes)   # get the 13 complete cases
# cc(nhanes$bmi) # extract complete bmi
}
