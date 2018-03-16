#' @include GarmaSim.R
NULL
#' \emph{GARMA} simulation fit class.
#'
#' An S4 class to represent a \emph{garma} simulation fitted object.
#'
#'This S4 class defines the basic structre of a \emph{GARMA}
#'simulation fitted object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@slot garma An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@slot allow.parallel  Logical \bold{TRUE/FALSE} indicating
#' whether parallel computation via the foreach package
#' should be used. The default value is \code{TRUE}. OBS:paralllel
#'  backend must be registered prior to calling \code{\link{GarmaSim}}.
#'
#'@slot seed Numeric, the seed to \code{set.seed()} for
#' replicable examples. Default value is \bold{123}.
#'
#'@slot errorhandling Character, either 'try' or 'pass'
#'
#'@slot n.try Positive integer. If \code{errorhandling = 'try'},
#'this specifies the number of attempts in the algorithm.
#'
#'@slot control List. This is passed to the garmaFit2 function. The options
#'are given by  \code{\link[gamlss.util]{garmaFit}}.
#'
#'@slot print.out A matrix with values to be used by \code{print}.
#'The column represents a simulation (\code{nmonte}) and the rows
#'the estimated parameters.
#'
#'@slot plot.out  A list with two elements: \describe{
#'   \item{db}{A data frame with columns: parameter,
#'   value and label (grouping variable used for some plots) 
#'   If \code{nmonte = 1} there is also variable column. Each row
#'   represents a statistic computed for a given block length, 
#'   \code{nmonte} and parameter.} 
#'   \item{db2}{A data frame with columns parameter, variable, 
#'   value and label and rows the model parameters. The values 
#'   are the simulated parameter values or parameter true-values.}
#' } with values to be used by \code{plot}.
#'
#'@slot value A list, where each element is a fitted
#'\emph{GARMA} object.
#'




setClass("GarmaFit",
         slots = list(garma = 'GarmaSim',
                      allow.parallel = 'logical', seed = 'numeric',
                      errorhandling = 'character' , n.try = 'numeric',
                      control = 'list',print.out = 'matrix',
                      plot.out = 'list', value = 'list'),
         prototype = list(allow.parallel = TRUE,
                          seed = 123, errorhandling = 'try',
                          n.try = 5L, control = list(iter.max=1000)),
         validity = function(object) {
           if(!object@errorhandling %in% c("try","pass"))
             return("unrecognized value of 'errorhandling',
		accepted values are: 'try' and 'pass'")

           if((object@errorhandling == "try")&&(object@n.try <= 0))
             return("'n.try' must be a positive integer")

           if(object@allow.parallel)
             if (foreach::getDoParRegistered()==FALSE)
               return("parallel backend must be registered")
           return(TRUE)
         }

)
