#' @include GarmaSim.R
NULL
#' \emph{GARMA} simulation bootstrap class.
#'
#' An S4 class to represent a \emph{garma} simulation bootstraped object.
#'
#'This S4 class defines the basic structre of a \emph{GARMA}
#'simulation bootstraped object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@slot sim An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#' @slot l \code{l} is the fixed block length used in generating the replicate
#'time series

#'
#'@slot R A positive integer giving the number of
#' bootstrap replicates required.
#'
#' @slot allow.parallel  Logical \bold{TRUE/FALSE} indicating
#' whether parallel computation via the foreach package
#' should be used. The default value is \code{TRUE}. OBS:paralllel
#'  backend must be registered prior to calling \code{\link{GarmaSimBoot}}.
#'
#'@slot seed Numeric, the seed to \code{set.seed()} for
#' replicable examples. Default value is \bold{123}.
#'
#'@slot errorhandling Character, either 'try' or 'pass'
#'
#'@slot n.try Positive integer. If \code{errorhandling = 'try'},
#'this specifies the number of attempts in the algorithm.
#'
#'@slot boot.function A function to summarise the bootstrap replicates.
#'The default function returns 0. Be aware that this is not a problem, as by
#'default the mean values is already being returned. This is useful if the user
#'wants to specify a quantity not being reported, as an example consider the
#'0.2 quantile.
#'
#'@slot control List. This is passed to the garmaFit2 function. The options
#'are given by  \code{\link[gamlss.util]{garmaFit}}.
#'
#'@slot print.out A data.frame with values to be used by \code{print}.
#' If sim@nmonte > 1 the rows display the the statitic computed 
#' for each parameter in the model for all \code{l} and \code{nmonte}.
#' The columns exhibit some chosen statistics that are 
#' computed with the resampled values, such as the original 
#' value (\code{l =n}), the estimated bias, the bias corrected
#' parameter estimate, the std. error etc. Though, if
#' sim@nmonte == 1, no statistic is computed and the estimated 
#' parameters in each block length is reported. 
#'
#'@slot plot.out  A list with two elements: \describe{
#'   \item{db}{A data frame with columns: length, parameter,
#'   variable (original, bias corrected etc...), value and 
#'   label (grouping variable used for some plots). Each row
#'   represents a statistic computed for a given block length, 
#'   \code{nmonte} and parameter. Note that if \code{nmonte = 1}
#'   no summarisation is performed and the original parameter
#'   values are returned.}
#'   \item{db2}{A data frame with columns parameter, variable, 
#'   value and label and rows the model parameters. The values 
#'   are the simulated parameter values or parameter true-values.}
#' } with values to be used by \code{plot}.
#'
#'@slot value A list, where each element is a fitted
#'\emph{GARMA} object.
#'
setClass("GarmaSimBoot",
         slots = list(sim = 'GarmaSim', l = 'numeric',
                      R = 'numeric', allow.parallel = 'logical',
                      seed = 'numeric', errorhandling = 'character',
                      n.try = 'numeric', boot.function = 'function',
                      control = 'list', print.out = 'data.frame',
                      plot.out = 'list', value = 'list'),
         prototype = list( R = 100L, l = integer(0),
                           allow.parallel = TRUE,seed = 123,
                           errorhandling = 'try', n.try = 5L,
                           boot.function = function(x) 0,
                           control = list(iter.max=1000)),
         validity = function(object) {
           if(identical(object@l,integer(0)))
             return("'l' must be a positive integer smaller
		than nsteps.")

           if((object@l <= 0)||(object@l > object@sim@nsteps))
             return("'l' must be a positive integer smaller
		than nsteps.")

           if(object@R <=0)
             return("'R' must be a positive integer.")

           if(!object@errorhandling %in% c("try","pass"))
             return("unrecognized value of 'errorhandling',
		accepted values are: 'try' and 'pass'.")

           if((object@errorhandling == "try")&&(object@n.try <= 0))
             return("'n.try' must be a positive integer.")

           if(object@allow.parallel)
             if (foreach::getDoParRegistered()==FALSE)
               return("parallel backend must be registered.")
           return(TRUE)
         }

)
