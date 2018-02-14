#' @include GarmaFit.R
NULL


#' GARMA fit  class generic function.
#'
#' A generic function working as a class constructor for the GarmaFit class.
#'
#' This generic function is designed to be a constructor for the
#' GarmaFit class.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@param garma An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@param \dots Further arguments that specify the GarmaFit object,
#'check the \bold{slots} arguments.
#'


setGeneric("GarmaFit", function(garma, ...) {
  standardGeneric("GarmaFit")
})



#'@describeIn GarmaFit
#'
#' GarmaFit Class constructor.
#'
#' A method to create an object of the GarmaFit class.
#'
#' This method creates an object of the GarmaFit class.
#'
#'@param garma An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@param \dots Further arguments that specify the GarmaFit object,
#'check the \bold{slots} arguments.



setMethod(f="GarmaFit",
          definition = function(garma, ...) {
            obj <- new("GarmaFit",garma = garma, ...)

            `%op%` <- if(obj@garma@allow.parallel)
              `%dorng%` else `%do%`
            set.seed(obj@seed)

            out <- foreach(i=obj@garma@value,
                           .packages='gamlss')%op%{
                             switch(obj@errorhandling
                                    ,"try"= {
                                      fit <- NULL
                                      attempt <- 1
                                      while( is.null(fit) && attempt <= obj@n.try ) {
                                        attempt <- attempt + 1
                                        try(fit <- dboot::garmaFit2(
                                          formula = formula(y ~. -1) ,
                                          order = obj@garma@order,
                                          data = as.data.frame(cbind(y = i$yt,
                                                                     tail(obj@garma@spec@X,
                                                                          obj@garma@nsteps) )),
                                          family = obj@garma@spec@family,
                                          control = obj@control
                                        )
                                        )
                                      }

                                    }
                                    ,"pass"={fit <- dboot::garmaFit2(
                                      formula = formula(y ~. -1) ,
                                      order = obj@garma@order,
                                      data = as.data.frame(cbind(y = i$yt,
                                                                 tail(obj@garma@spec@X,
                                                                      obj@garma@nsteps))),
                                      family = obj@garma@spec@family,
                                      control = obj@control)

                                    }
                             )
                             fit
                           }

            slot(obj,"value") <- out

            print.out <- if(obj@garma@nmonte > 1) {
              if(all(obj@garma@order == c(0,0)))
              {db <- sapply(obj@value, function(j) j$mu.coef)
              rownames(db) <- paste("beta.",rownames(db),sep="")
              db
              } else
              {
                db <- sapply(obj@value, function(j) j$coef)
                ind <- match("beta.null.vector",dimnames(db)[[1]])
                db <- if(!is.na(ind)) db[-ind,] else db
                ind <- match("beta.null.vector",dimnames(db)[[2]])
                db <- if(!is.na(ind)) db[,-ind] else db
                db
              }
            } else
               obj@value[[1]]$coef

            slot(obj,"print.out") <- as.matrix(print.out)

            plot.out <- {

              create.labels <- function(data,labels =
                                          c("beta","phi","theta","sigma")) {
                data$label <- 0
                for (i in labels) {
                  data$label[grep(i,data$parameters)] <- i
                }
                return(data)
              }
              prep.db <- function(data) {
                data <- reshape2::melt(data)
                data <- data[,-2]
                colnames(data) <- c("parameters","value")
                data <- create.labels(data)
                return(data)
              }


              create.parameters <- function(x) {
                par1 <- if(identical(grep("beta", rownames(x@print.out)),
                                     integer(0))) NULL else x@garma@spec@beta.x
                par2 <- if(identical(grep("phi", rownames(x@print.out)),
                                     integer(0))) NULL else x@garma@spec@phi
                par3 <- if(identical(grep("theta", rownames(x@print.out)),
                                     integer(0))) NULL else x@garma@spec@theta
                par4 <- if(identical(grep("sigma", rownames(x@print.out)),
                                     integer(0))) NULL else sqrt(x@garma@spec@sigma2)
                pars <- c(par1,par2,par3,par4)
                return(pars)
              }

              db <- prep.db(obj@print.out)

              db2 <- data.frame(parameters = levels(db$parameters),
                                variable = factor("true.value"),
                                value = create.parameters(obj))
              db2 <- create.labels(db2)

              list(db = db, db2=db2)

            }

            slot(obj,"plot.out") <- plot.out

            return(obj)

          })

#'@describeIn GarmaFit
#'
#' Print method for a GarmaFit object.
#'
#'@param x An object of the GarmaFit class, as provided
#'by \code{\link{GarmaFit}}.
#'

setMethod(f="print",
          signature = "GarmaFit",
          definition = function(x) {
            cat(" -------------------------------------------------------\n",
                "A",paste(x@garma@spec@family,"-Garma(",
                          x@garma@order[1],",",x@garma@order[2],")",
                          sep=""),"simulation fitted object: \n\n",

                if(x@garma@order[1] != 0) paste(
                  if(x@garma@order[1] == 1) "phi" else
                    paste("phi",1:x@garma@order[1],
                          " = ",sep=""), x@garma@spec@phi,"\n"),

                if(x@garma@order[2] != 0) paste(
                  if(x@garma@order[2] == 1) "theta" else
                    paste("theta",1:x@garma@order[2],
                          " = ",sep=""), x@garma@spec@theta,"\n"),

                if(TRUE%in%(colnames(x@garma@spec@X) != "null.vector"))
                  paste(paste("beta.",colnames(x@garma@spec@X), " = ",
                              sep=""), x@garma@spec@beta.x,"\n"),

                if(TRUE%in%(x@garma@spec@mu0 != 0)) paste(paste("mu0[", #altersr
                                                                1:length(x@garma@spec@mu0),"] = ",sep=""),
                                                          x@garma@spec@mu0,"\n"),

                if((x@garma@spec@family == "GA")&&(TRUE%in%(x@garma@spec@sigma2 != 0)))
                  paste("sigma2 = ", paste(x@garma@spec@sigma2),"\n"),

                "Number of Monte Carlo Simulations ('nmonte') = ",x@garma@nmonte,"\n",
                "Time Series Length ('nsteps') = ",x@garma@nsteps,"\n",
                "Burn-in ('burnin') = ",x@garma@burnin,"\n",
                "parallel = ",x@allow.parallel,"\n",
                "errorhandling = ",x@errorhandling,"\n",
                "n.try = ",x@n.try,"\n",

                "-------------------------------------------------------\n",
                "\n",
                "Estimated parameters: \n\n")

            x@print.out

          }
)


#'@describeIn GarmaFit
#'
#' Plot method for a GarmaFit object.
#'
#'@param x An object of the GarmaFit class, as provided
#'by \code{\link{GarmaFit}}.
#'
#'@param scales A character specifying if the scales should be free or not.
#'The default value is \code{"free"} and accepted values are:
#' \code{'fixed', 'free_x', 'free_y' and 'free'}.

setMethod(f="plot",
          signature = "GarmaFit",
          definition = function(x,scales = "free",...) {
            if(!scales%in%c('fixed', 'free_x', 'free_y', 'free'))
              stop("invalid value of scales, accepted values:
                   'fixed', 'free_x', 'free_y', 'free' ")

            db <- x@plot.out[['db']]
            db2 <- x@plot.out[['db2']]

            if((x@garma@nmonte > 1))
            {
              names(db2)[1] <- "true.value"

              if((length(x@garma@spec@beta.x) > 1)||
                 (x@garma@order[1] >1)||(x@garma@order[2] >1)) {
                ggplot(db, aes(value, fill = parameters)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  facet_grid(label~.,scales=scales) +
                  ggtitle(paste("Simulated", paste(x@garma@spec@family,"-Garma(",
                                                   x@garma@order[1],",",x@garma@order[2],"): ",	sep=""),
                                " Coefficients Distribution")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5))

              } else {
                ggplot(db, aes(value, fill = parameters)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  ggtitle(paste("Simulated", paste(x@garma@spec@family,"-Garma(",
                                                   x@garma@order[1],",",x@garma@order[2],"): ",	sep=""),
                                " Coefficients Distribution")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5))

              }

            } else {

              db$variable <- factor("estimated.value")
              db <- rbind(db,db2)
              db$text <- paste(round(db$value,2))
              ggplot(db, aes(x = parameters, y = value, colour = parameters,
                             shape= variable))	+ geom_point(size=10) +
                geom_text(aes(label=text),size = 4 , vjust=0.4) +
                scale_shape_manual(values=c(0,1), name="type.param") +
                labs(colour = "parameter",x = "parameter") +
                ggtitle(paste(paste(x@garma@spec@family,"-Garma(",
                                    x@garma@order[1],",",x@garma@order[2],"): ",	sep=""),
                              " Estimated vs. True Value parameters")) +
                theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                              face="bold.italic", colour="gray26",hjust=0.5))

            }

          }
)

#'@describeIn GarmaFit
#'
#' Summary method for a GarmaFit object.
#'
#'@param object An object of the GarmaFit class, as provided
#'by \code{\link{GarmaFit}}.
#'


setMethod(f="summary",
          signature = "GarmaFit",
          definition = function(object,...) {

            cat(" -------------------------------------------------------\n",
                "A",paste(object@garma@spec@family,"-Garma(",
                          object@garma@order[1],",",object@garma@order[2],") ",
                          sep=""),"simulation object: \n\n",
                "Number of Monte Carlo Simulations ('nmonte') = ",object@garma@nmonte,"\n",
                "Time Series Length ('nsteps') = ",object@garma@nsteps,"\n",
                "Burn-in ('burnin') = ",object@garma@burnin,"\n",
                "parallel = ",object@garma@allow.parallel,"\n",
                "-------------------------------------------------------\n")

            if(object@garma@nmonte > 1)
            {
              db <- plyr::ddply(object@plot.out$db,
                                "parameters",'summarise',Min. = min(value),
                                "1st Qu." = quantile(value,(.25),type = 1),
                                "Median." = quantile(value,(.50),type = 1),
                                "Mean." = mean(value),
                                "3rd Qu." = quantile(value,(.75),type = 1),
                                "Max." = max(value))
              db <-cbind(db,true.value=object@plot.out$db2$value)
              return(db)
            } else
              db <- object@plot.out$db2[c('parameters','value')]
            names(db)[2] <- 'true.value'
            db <- cbind(db,estimate=object@plot.out$db$value)
            return(db)
          }
)
