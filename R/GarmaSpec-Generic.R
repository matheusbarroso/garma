#' @include GarmaSpec.R
NULL

#' GARMA specification class generic function.
#'
#' A generic function working as a class constructor for the GarmaSpec class.
#'
#' This generic function is designed to be a constructor for the
#' GarmaSpec family.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@param family A character vector specifying the
#'family of the specification object. Accepted values
#'are: "Po" for poisson and "GA"for gamma families.
#'
#'@param \dots Further arguments that specify the model,
#'check the \bold{slots} arguments.





setGeneric("GarmaSpec", function(family, ...) {
  standardGeneric("GarmaSpec")
})




#'@describeIn GarmaSpec
#'
#' GarmaSpec Class constructor.
#'
#' A method to create an object of the GarmaSpec class.
#'
#' This method creates an object of the GarmaSpec class.
#'
#'@param family A character vector specifying the
#'family of the specification object. Accepted values
#'are: "Po" for poisson and "GA"for gamma families.
#'
#'@param \dots Further arguments that specify the model.
#'
#'@examples #PoissonSpec 
#'GarmaSpec("PO",phi=c(1,1),mu0=c(1,1))

#'#GammaSpec
#'
#'GarmaSpec("GA",phi=c(1,1),mu0=c(1,1))



setMethod(f="GarmaSpec",

          definition = function(family, ...) {
          if(!family %in% c("PO","GA"))
            stop("Invalid value of family, accepted values 
                 are: 'PO' and 'GA'")
            switch(family,"PO" = {
              new("PoissonSpec",family = family, ...)

            }
            , "GA" = {
              new("GammaSpec",family = family, ...)

            }
            )
          }

)
