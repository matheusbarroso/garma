#' @include GarmaSpec.R
NULL
#' \emph{Poisson}-GARMA specification class.
#'
#' An S4 class to represent a \emph{Poisson} garma specification object.
#'
#'This is S4 class defines the basic structre of a \emph{Poisson} GARMA
#'speficication object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@slot family A character with the tag "PO".
#'
#'@slot alpha Numeric, specifying the \bold{offset} parameter.
#'
#'@slot mu0 A numeric vector with length equal to the \bold{max.order}
#'of the \emph{Poisson}-GARMA model. See \code{\link{GarmaSpec}} for
#'further details.
#'

setClass("PoissonSpec",
         slots = list(alpha = "numeric",
                      mu0 = "numeric"),
         prototype = list(family = "PO", alpha = 0.1,
                          mu0 = 1, y0 = 1),
         validity = function(object) {
           if(object@family != "PO")
             return("The PoissonSpec class is only valid
			for family = 'PO'.")

           if(object@alpha < 0)
             return("The offset term must be positive.")

           if(TRUE%in%(object@mu0 < 0))
             return("The mu0 term must be positive.")

           if(length(object@mu0) < max(length(object@phi),
                                       length(object@theta)))
             return("The length of mu0 should be
		of the highest order of the Garma model")

           return(TRUE)
         },
         contains = "GarmaSpec"
)
