#' @include GarmaSimBoot.R
NULL

#' \emph{GARMA} MBB simulation confidence interval class.
#'
#' An S4 class to represent a \emph{garma} MBB simulation 
#' confidence interval object.
#'
#'This S4 class defines the basic structre of a \emph{GARMA}
#'MBB simulation confidence interval object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@slot garma.boot An object of the \bold{GarmaSimBoot}
#' class, as provided by \code{\link{GarmaSimBoot}}.
#'
#'@slot allow.parallel  Logical \bold{TRUE/FALSE} indicating
#' whether parallel computation via the plyr (foreach) package
#' should be used. The default value is \code{TRUE}. OBS:paralllel
#'  backend must be registered prior to calling
#'   \code{\link{ConfidenceInterval}}.
#'
#'@slot conf A scalar containing the confidence level 
#'of the required interval(s).
#'
#'@slot summary.out A data.frame with the mean 
#'(w.r.t. \code{nmonte}) values for all block lengths,
#'parameters and ci's.
#' 
#'@slot plot.out A \code{data.frame} with the computed 
#'ci's for each parameter, block length and \code{nmonte}.
#'
#'@slot coverage.out A \code{data.frame} with the computed 
#'ci's coverage rates for each block length , parameter 
#'and type of ci (Normal, Bias Corrected, Basic, Percentile)
#'
#'@slot value A list, where each element is a data.frame
#'with the computed ci's for each parameter in each Monte 
#'Carlo simulation through the MBB resampling process. 
#'The number of elements in the list is equivalent to 
#'\code{l} in \code{\link{GarmaSimBoot}}.
#'




setClass("ConfidenceInterval",
         slots = list(garma.boot = 'GarmaSimBoot',
                      allow.parallel = 'logical',
                      conf = 'numeric',
                      summary.out = 'data.frame',
                      plot.out = 'data.frame',
                      coverage.out = 'data.frame',
                      value = 'list'),
         prototype = list(allow.parallel = TRUE,
                          conf =  0.95
                          ),
         validity = function(object) {
           if(object@allow.parallel)
             if (foreach::getDoParRegistered()==FALSE)
               return("parallel backend must be registered")
           return(TRUE)
           
           if(!(object@conf >= 0)&&(object@conf <= 1))
             return("'conf' should be larger than or equal 
                    to zero and smaller than or euqual to one.
                    That is: 0 <= conf <= 1")
         }
         
)
