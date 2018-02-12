#' @include GarmaSpec.R
NULL
#' \emph{Gamma}-GARMA specification class.
#'
#' An S4 class to represent a \emph{Gamma} garma specification object.
#'
#'This is S4 class defines the basic structre of a \emph{Gamma} GARMA
#'speficication object.
#'
#'@slot family A character with the tag "GA".
#'
#'@slot sigma2 Numeric, specifying the \bold{sigma2} parameter.
#'
#'@slot mu0 A numeric vector with length equal to the \bold{max.order}
#'of the \emph{Gamma}-GARMA model. See \code{\link{GarmaSpec}} for
#'further details.
#'
#'
setClass("GammaSpec",
         slots = list(sigma2 =  "numeric",
                      mu0 = "numeric"),
         prototype = list(family = "GA", sigma2 = 1,
                          mu0 = 10, y0 = 1),
         validity = function(object) {
           if(object@family != "GA")
             return("The PoissonSpec class is only valid
                    for family = 'GA'.")

           if(object@sigma2 < 0)
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
