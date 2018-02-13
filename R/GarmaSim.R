#'@include GarmaSpec.R
NULL

#' \emph{GARMA} simulation class.
#'
#' An S4 class to represent a \emph{garma} simulation object.
#'
#'This S4 class defines the basic structre of a \emph{GARMA}
#'simulation object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@slot spec An object of the \bold{GarmaSpec} class, as provided
#'by \code{\link{GarmaSpec}}.
#'
#' @slot nmonte A positive integer, specifying the number of Monte Carlo
#' simulations to perform, the default value is \bold{1000}.
#'
#' @slot nsteps A numeric vector with the number of steps in the
#' Garma model simulation, that is, the length of the time series
#' to simulate. The default value is \bold{100}.
#'
#' @slot burnin A numeric vector indicating the number of burn in
#' observations. If you want to generate only \code{nsteps} this
#' argument should be set equal to zero. Otherwise, provide a positive
#' integer. The default value is \bold{1000}.
#'
#' @slot allow.parallel  Logical \bold{TRUE/FALSE} indicating
#' whether parallel computation via the foreach package
#' should be used. The default value is \code{TRUE}. OBS:paralllel
#'  backend must be registered prior to calling \code{\link{GarmaSim}}.
#'
#'@slot seed Numeric, the seed to \code{set.seed()} for
 #' replicable examples. Default value is \bold{123}.


setClass("GarmaSim",
         slots = list(spec = "GarmaSpec", nmonte = "numeric",
                      nsteps = "numeric", burnin = "numeric",
                      allow.parallel = 'logical', seed = "numeric",
                      print.out ="matrix", plot.out ="data.frame",
                      summary.out = "list",
                      value = "list", order = "numeric"),
         prototype = list(nmonte = 1000, nsteps = 100,
                          burnin = 1000, allow.parallel = TRUE,
                          seed = 123),
         validity = function(object) {
           if(!validObject(object@spec))
             return("spec must be a valid GarmaSpec object")

           if(object@nmonte < 1)
             return("The number of Monte Carlo simulations
		must be a positive integer")

           if(object@nsteps < 1)
             return("The number of steps (length of the
		simulated series) must be a positive integer")

           if(object@burnin < 0)
             return("The burnin must be zero or a positive integer")

           return(TRUE)
         }
)
