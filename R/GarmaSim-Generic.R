#' @include GarmaSim.R
NULL


#' GARMA simulation class generic function.
#'
#' A generic function working as a class constructor for the GarmaSim class.
#'
#' This generic function is designed to be a constructor for the
#' GarmaSim class.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@param spec An object of the GarmaSpec class, as provided
#'by \code{\link{GarmaSpec}}.
#'
#'@param \dots Further arguments that specify the GarmaSim object,
#'check the \bold{slots} arguments.



setGeneric("GarmaSim", function(spec, ...) {
  standardGeneric("GarmaSim")
})

#'@describeIn GarmaSim
#'
#' GarmaSim Class constructor.
#'
#' A method to create an object of the GarmaSim class.
#'
#' This method creates an object of the GarmaSim class.
#'
#'@param spec An object of the GarmaSpec class, as provided
#'by \code{\link{GarmaSpec}}.
#'
#'@param \dots Further arguments that specify the simulation object.
#'Check the \bold{Slots} for further details.
#'
#'@examples ##Some Specs/Simulations for different outputs:
#'
#'\dontrun{
#'#Incompatible dimensions x and beta.x:
#'ex1 <- GarmaSim(
#'GarmaSpec("GA", 
#'phi = c(0.2),
#'theta = c(1,2,3),
#'mu0 = 1:3,
#'X = as.matrix(
#'data.frame(
#'intercept = rep(4,1100),
#'x1 = rep(3,1100)))),
#'nmonte = 2)
#'
#'#Incompatible number of rows
#'ex2 <- GarmaSim(
#'GarmaSpec("GA", 
#'phi = c(0.2),
#'theta = c(1,2,3),
#'mu0 = 1:3,
#'beta.x = c(1,2),
#'X = as.matrix(
#'data.frame(
#'intercept = rep(4,1100),
#'x1 = rep(3,1100)))),
#'nmonte = 2)
#'}
#'
#'ex3 <- GarmaSim(
#'GarmaSpec("PO",
#'beta.x = c(0.1,1),
#'X = as.matrix(
#'data.frame(
#'intercept = rep(10,1100),
#'x1 = c(rep(7,100),rep(2,1000))))),
#'nmonte = 1, allow.parallel = TRUE) 
#'
#'ex4 <- GarmaSim(
#'GarmaSpec("PO", 
#'phi = c(0.2),
#'theta = c(1,2,3),
#'mu0 = 1:3),
#'nmonte = 2)
#'
#'ex5 <- GarmaSim(
#'GarmaSpec("PO", 
#'mu0 = 10,
#'phi = 0.15, 
#'X = as.matrix(
#'data.frame(
#'  x1 = rep(10,101)))),
#'burnin = 0,
#'nmonte = 10)
#'
#'ex6 <- GarmaSim(
#'GarmaSpec("GA", 
#'phi = c(0.2),
#'theta = c(1,2,3),
#'mu0 = 1:3,
#'beta.x = c(1,1),
#'X = as.matrix(
#'data.frame(
#'intercept = rep(4,1103),
#'  x1 = rep(3,1103)))),
#'nmonte = 2)
#'
#'ex7 <- GarmaSim(
#'GarmaSpec("GA", 
#'phi = c(0.2),
#'theta = c(1,2,3),
#'mu0 = 1:3,
#'beta.x = c(1,2),
#'X = as.matrix(
#'data.frame(
#'intercept = rep(4,1103),
#' x1  =rep(3,1103)))),
#'nmonte = 10)
#'
#'ex8 <- GarmaSim(
#'GarmaSpec("PO", 
#'phi = c(0.5, 0.15),
#'mu0 = c(1000, 1100),
#'beta.x = c(1, 1),
#'X = as.matrix(
#'data.frame(
#'intercept = rep(1, 1102),
#'x1 = c(rep(7, 100),
#'rep(2,1002))))),
#'nmonte = 10))
#'
#'ex9 <- GarmaSim(
#'GarmaSpec("GA"), 
#'nmonte = 10)
#'
#'ex10 <- GarmaSim(
#'GarmaSpec("PO"), 
#'nmonte = 10)
#'
#'
#'
#' ##Example of the GarmaSim methods:
#'
#'print(ex1)
#'plot(ex1)
#'summary(ex1)

setMethod(f="GarmaSim",
          definition = function(spec, ...) {
            obj <- new("GarmaSim",spec = spec, ...)
            slot(obj,"order") <- c(if(length(spec@phi) > 1)
              length(spec@phi) else if(spec@phi==0) 0 else 1,
              if(length(spec@theta) > 1)
                length(spec@theta) else if(spec@theta==0) 0 else 1)

            if(obj@allow.parallel)
              if (foreach::getDoParRegistered()==FALSE)
                stop("parallel backend must be registered or
                     set allow.parallel = FALSE")

            `%op%` <- if(obj@allow.parallel) `%dorng%` else `%do%`

            total.length <- obj@nsteps + obj@burnin

            if(all(dim(obj@spec@X) == 0))
              slot(slot(obj,"spec"),"X") <- as.matrix(data.frame(
                null.vector = numeric(total.length +
                                        max(obj@order)
                )))

            else {
              if(dim(obj@spec@X)[2] != length(obj@spec@beta.x))
                stop("X and beta.x dimensions are incompatible")
              if(dim(obj@spec@X)[1] != total.length+max(obj@order))
                stop("X must have 'nsteps'+'burnin'+'max. order' rows")
            }

            B0 <- if(is.null(names(obj@spec@beta.x))|is.null(
              colnames(obj@spec@X)))
              obj@spec@X%*%obj@spec@beta.x
            else {
              if(!all(names(obj@spec@beta.x) %in% colnames(obj@spec@X)))
                stop("The names in 'beta.x' do not match the ones in 'X',
                     please provide matching arguments or an unnamed vector.
                     In the latter case no matching will be done.") else
                       obj@spec@X[,match(names(obj@spec@beta.x),
                                         colnames(obj@spec@X))]%*%obj@spec@beta.x
            }


            set.seed(obj@seed) # maybe unecessary, just check if there is a seed or set one random.
            link.function <- function(x) switch(spec@family,
                                                ,"PO" = {max(log(x),obj@spec@alpha) }
                                                ,"GA" = {log(x)}
            )

            out <- foreach(n=seq_len(obj@nmonte))%op% { #obj@nmonte

              y0 <- switch(obj@spec@family
                           ,"PO" = { rpois(length(obj@spec@mu0),length(obj@spec@mu0))}
                           ,"GA" = {rgamma(length(obj@spec@mu0),shape=1/obj@spec@sigma2,
                                           rate=1/(obj@spec@sigma2*length(obj@spec@mu0)))}
              )

              db0 <- if(all(obj@order == c(0,0))) data.frame() else
                data.frame(indext = seq.int(to = 0,
                                            length.out = length(obj@spec@mu0)),
                           mu.t=obj@spec@mu0,yt=y0)
              shift <- nrow(db0)
              db <- data.frame(indext = seq_len(total.length),
                               mu.t=rep(0,total.length),yt=rep(0,total.length))
              db <- rbind(db0,db)
              if(!all(dim(obj@spec@X) == 0)) db <- cbind(db,B0)


              for (k in seq_len(total.length)) {

                actual.position <- k + shift

                db$mu.t[actual.position] <- exp(B0[actual.position] +
                                                  if(obj@order[1]==0) 0 else {
                                                    sum(sapply(seq_len(obj@order[1]), function(i)
                                                      obj@spec@phi[i]*(link.function(db$yt[actual.position-i])-
                                                                         B0[actual.position-i])))
                                                  } +
                                                  if(obj@order[2]==0) 0 else {
                                                    sum(sapply(seq_len(obj@order[2]), function(i)
                                                      obj@spec@theta[i]*(link.function(db$yt[actual.position-i]) -
                                                                           link.function(db$mu.t[actual.position-i]))))
                                                  }
                )

                db$yt[actual.position] <- switch(spec@family,
                                                 ,"PO" = { rpois(1,db$mu.t[actual.position])}
                                                 ,"GA" = {rgamma(1,shape=1/obj@spec@sigma2,
                                                                 rate=1/(obj@spec@sigma2*db$mu.t[actual.position]))}
                )


              }
              tail(db,obj@nsteps)
             # db[seq.int(from = nrow(db)-obj@nsteps + 1,
                    #     to = nrow(db)),]

            }
            slot(obj,"value") <- out


            print.out <- {
              db <- sapply(obj@value, function(j) j$yt)
              colnames(db) <- paste('yt',seq_len(obj@nmonte),sep="_")
              rownames(db) <- seq_len(obj@nsteps)
              db
            }

            slot(obj,"print.out") <- print.out


            plot.out <- {
              db <- sapply(obj@value, function(j) j$yt)
              colnames(db) <- paste('yt',seq_len(obj@nmonte),sep="_")
              rownames(db) <- seq_len(obj@nsteps)

              as.data.frame(db)
            }

            slot(obj,"plot.out") <- plot.out

            summary.out <- {
              db <- sapply(obj@value, function(j) c(min(j$yt),
                                                       mean(j$yt),
                                                       max(j$yt)))

              db1 <- as.data.frame(t(apply(db,1,quantile,
                                           probs=c((1-0.95)/2,0.5,0.95+(1-0.95)/2),
                                           type=1)))

              db1 <- round(db1,4)

              db <- round(apply(db,1,mean),4)

              rownames(db1) <- NULL

              list(db=db,db1=db1)
            }
            slot(obj,"summary.out") <- summary.out

            return(obj)
          }
)

#'@describeIn GarmaSim
#'
#' Print method for a GarmaSim object.
#'
#'@param x An object of the GarmaSim class, as provided
#'by \code{\link{GarmaSim}}.
#'


setMethod(f="print",
          signature = "GarmaSim",
          definition = function(x) {
            cat(" -------------------------------------------------------\n",
                "A",paste(x@spec@family,"-Garma(",
                          x@order[1],",",x@order[2],")",
                          sep=""),"simulation object: \n\n",
                if(x@order[1] != 0) paste(
                  if(x@order[1] == 1) "phi" else
                    paste("phi",1:x@order[1],
                          " = ",sep=""), x@spec@phi,"\n"),

                if(x@order[2] != 0) paste(
                  if(x@order[2] == 1) "theta" else
                    paste("theta",1:x@order[2],
                          " = ",sep=""), x@spec@theta,"\n"),

                if(TRUE%in%(colnames(x@spec@X) != "null.vector"))
                  paste(paste("beta.",colnames(x@spec@X), " = ",
                              sep=""), x@spec@beta.x,"\n"),

                if(TRUE%in%(x@spec@mu0 != 0)) paste(paste("mu0[", #altersr
                                                          1:length(x@spec@mu0),"] = ",sep=""),
                                                    x@spec@mu0,"\n"),

                if((x@spec@family == "GA")&&(TRUE%in%(x@spec@sigma2 != 0)))
                  paste("sigma2 = ", paste(x@spec@sigma2),"\n"),
                "Number of Monte Carlo Simulations ('nmonte') = ",x@nmonte,"\n",
                "Time Series Length ('nsteps') = ",x@nsteps,"\n",
                "Burn-in ('burnin') = ",x@burnin,"\n",
                "parallel = ",x@allow.parallel,"\n",
                "-------------------------------------------------------\n",
                "\n")

            if(x@nmonte > 1)
              x@print.out else
                x@value[[1]]$yt

          }
)


#'@describeIn GarmaSim
#'
#' Print method for a GarmaSim object default plot.
#'
#'
#'@param x An object of the GarmaSim class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@param confInt A numeric value, \eqn{0 < confInt <1},
#'giving the desired confidente interval size for the plot.
#'

setMethod(f="plot",
          signature = "GarmaSim",
          definition = function(x,confInt=0.95,...) {

            if(x@nmonte > 1)
            {db <- x@plot.out
            db <- cbind(data.frame(index_t = seq_len(x@nsteps),
                                   mean_yt =
                                     apply(db,1,mean)),
                        as.data.frame(t(apply(db,1,quantile,
                                              probs=c((1-confInt)/2,0.5,
                                                      confInt+(1-confInt)/2),
                                              type=1))))

            colnames(db) <- c("index_t","mean_yt","lim.inf","median","lim.sup")
            ggplot(db, aes(x = index_t, y = mean_yt)) +
              geom_ribbon(aes(ymin = lim.inf, ymax = lim.sup,
                              alpha = 0.2)) +  #geom_line(aes(y=median), colour="orange") +
              theme_bw() + geom_line() +
              ggtitle(paste("Simulated", paste(x@spec@family,"-Garma(",
                                               x@order[1],",",x@order[2],"): ",	sep=""),"Mean values and",
                            paste(round(confInt*100,4),"%",sep="")," confidence
		interval through time")) +
              theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                            face="bold.italic", colour="gray26",hjust=0.5),
                    legend.position="none")

            } else
              ggplot(data.frame(index_t = seq_len(x@nsteps), yt= x@value[[1]]$yt),
                     aes(x = index_t, y = yt)) + theme_bw() +geom_line() +
              ggtitle(paste("Simulated", paste(x@spec@family,"-Garma(",
                                               x@order[1],",",x@order[2],")",	sep=""),"values through time")) +
              theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                            face="bold.italic", colour="gray26",hjust=0.5),
                    legend.position="none")

          }
)


#'@describeIn GarmaSim
#'
#' Summary method for a GarmaSim object default plot.
#'
#'@param object An object of the GarmaSim class, as provided
#'by \code{\link{GarmaSim}}.
#'



setMethod(f="summary",
          signature = "GarmaSim",
          definition = function(object,...) {

            cat(" -------------------------------------------------------\n",
                "A",paste(object@spec@family,"-Garma(",
                          object@order[1],",",object@order[2],") ",
                          sep=""),"simulation object: \n\n",
                "Number of Monte Carlo Simulations ('nmonte') = ",object@nmonte,"\n",
                "Time Series Length ('nsteps') = ",object@nsteps,"\n",
                "Burn-in ('burnin') = ",object@burnin,"\n",
                "parallel = ",object@allow.parallel,"\n",
                "-------------------------------------------------------\n")

            if(object@nmonte > 1)
            {db <- object@summary.out$db
            db1 <- object@summary.out$db1

            cat("Mean Monte Carlo mean values = ",db[2],"with distribution: \n")
            print(db1[2,])
            cat(" \n")
            cat("Mean Monte Carlo min values = ",db[1],"with distribution: \n")
            print(db1[1,])
            cat(" \n")
            cat("Mean Monte Carlo max values = ",db[3],"with distribution: \n")
            print(db1[3,])
            cat(" \n")
            db <- cbind(mean = db,db1)
            rownames(db) <- c("min","mean","max")
            db

            } else
              summary(object@value[[1]]$yt)

          }
)

