#' GARMA specification class.
#'
#' An S4 class to represent a garma specification object.
#'
#'This S4 class defines the basic structre of a GARMA
#'speficication object.
#'
#'@slot family A character vector specifying the
#'family of the specification object. Accepted values
#'are: "Po" for poisson and "GA"for gamma families.
#'
#' @slot beta.x A numeric vector with length of
#' the desired specification.
#' @slot phi A numeric vector with length of
#' the desired autoregressive term order specification.
#'
#' @slot theta A numeric vector with length of
#' the desired moving-average term order specification.
#'
#' @slot X A n x m matrix, where \code{n = nsteps + burnin + max.order}
#' and  \code{m = length(\strong{beta.x})}. Where \strong{nsteps}
#' is the number of simulations desired, \strong{burnin}
#' the number of burnin observations and \strong{max.order} is
#' equal to the highest order of the GARMA model (i.e.
#' \code{max(length(phi),length(theta)}.)
#'

setClass("GarmaSpec",
         slots = list( family = "character",beta.x = "numeric",
                       phi = "numeric", theta = "numeric",
                       X = "matrix"),
         prototype = list(beta.x = 1L,
                          phi = 0L, theta = 0L),
         validity = function(object) {

           return(TRUE)				}

)

