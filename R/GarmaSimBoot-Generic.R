#' @include GarmaSimBoot.R
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
#'@param sim  An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@param \dots Further arguments that specify the GarmaSimBoot object,
#'check the \bold{slots} arguments.
# #'@describeIn GarmaSimBoot


setGeneric("GarmaSimBoot", function(sim, ...) {
  standardGeneric("GarmaSimBoot")
})


#'@describeIn GarmaSimBoot
#'
#' GarmaSimBoot Class constructor.
#'
#' A method to create an object of the GarmaSimBoot class.
#'
#' This method creates an object of the GarmaSimBoot class.
#'
#'@param sim  An object of the \bold{GarmaSim} class, as provided
#'by \code{\link{GarmaSim}}.
#'
#'@param \dots Further arguments that specify the GarmaSimBoot object,
#'check the \bold{slots} arguments.



setMethod(f="GarmaSimBoot",
          definition = function(sim, ...) {
            obj <- new("GarmaSimBoot",sim = sim, ...)

            bootf <- function(data,order,family,control) {

              fit <- dboot::garmaFit2(
                formula = formula(y ~. -1) ,
                order = order,
                data = data,
                family = family,
                control = control
              )
              coeff <- fit$coef
              if(all(order == c(0,0)))
                coeff <- fit$mu.coef

              return(coeff)
            }

            h <- function(w) if( any( grepl( "<anonymous>:...", w) ) )
              invokeRestart( "muffleWarning" )

            out <- 	{
              lst <- lapply(obj@l, function(l){
                lapply(obj@sim@value, function(j) {
                  switch(obj@errorhandling
                         ,"try"= {
                           MBB <- NULL
                           attempt <- 1
                           while( is.null(MBB) && attempt <= obj@n.try ) {
                             attempt <- attempt + 1
                             try(MBB <- withCallingHandlers(
                               dboot::tsboot2(tseries =
                                                as.data.frame(cbind(y = j$yt,
                                                                    tail(obj@sim@spec@X,
                                                                         obj@sim@nsteps))),
                                              statistic = bootf,
                                              R = obj@R, l = l,
                                              allow.parallel = obj@allow.parallel,
                                              seed = obj@seed + attempt -1,
                                              packages = c("gamlss","dboot"),
                                              export = c("garmaFit2"),
                                              order = obj@sim@order,
                                              family = obj@sim@spec@family,
                                              control = obj@control),
                               warning = h )
                             )
                           }


                         }
                         ,"pass"={MBB <-  withCallingHandlers(
                           dboot::tsboot2(tseries =
                                            as.data.frame(cbind(y = j$yt,
                                                                tail(obj@sim@spec@X,
                                                                     obj@sim@nsteps) )),
                                          statistic = bootf,
                                          R = obj@R, l = l,
                                          allow.parallel = obj@allow.parallel,
                                          seed = obj@seed,
                                          packages = c("gamlss","dboot"),
                                          export = c("garmaFit2"),
                                          order = obj@sim@order,
                                          family = obj@sim@spec@family,
                                          control = obj@control),
                           warning = h )
                         }

                  )

                  MBB
                })
              })

              names(lst) <- paste("l=",obj@l,sep="")
              lst
            }

            slot(obj,"value") <- out

            clean <- plyr::llply(obj@value,function(j) {
              boot.function <- obj@boot.function
              plyr::llply(j,.parallel = obj@allow.parallel,
                          .fun = function(x) {
                            with(x,{
                              names(t) <- names(t0)
                              cbind(original = t0,
                                    bias = apply(t, 2L, mean, na.rm = TRUE) -  t0,
                                    "bias corrected" = 2*t0-(apply(t, 2L, mean, na.rm = TRUE)),
                                    " std. error" = apply(t,2L, sd,na.rm=TRUE),
                                    "Min." =  apply(t,2L, min,na.rm=TRUE),
                                    "1st Qu." = apply(t,2L, quantile,na.rm=TRUE, type =1, .25),
                                    "Median." = apply(t,2L, quantile,na.rm=TRUE, type =1, .50),
                                    "Mean." = apply(t, 2L, mean, na.rm = TRUE),
                                    "3rd Qu." = apply(t,2L, quantile,na.rm=TRUE, type =1, .75),
                                    "Max." =  apply(t,2L, max,na.rm=TRUE),
                                    boot.func = plyr::adply(.data = t,.margins = 2L,
                                                            .fun = boot.function, .id = NULL)
                              )
                            })
                          })
            }
            )

            print.out <- {

              out <- plyr::ldply(clean,function(j) {plyr::ldply(j,
                                                                function(i){
                                                                  i$parameter <- rownames(i)
                                                                  i })})
              names(out)[1] <- 'length'
              out[c(1,ncol(out),seq.int(2,ncol(out)-1))]

            }
            slot(obj,"print.out") <- print.out


            plot.out  <- plyr::ldply(clean,function(j) {plyr::ldply(j,
                                                                    function(i){
                                                                      i$parameter <- rownames(i)
                                                                      reshape2::melt(i,id.vars="parameter")
                                                                    })})

            slot(obj,"plot.out") <- plot.out


            return(obj)

          }
)

#'@describeIn GarmaSimBoot
#'
#' Print method for a GarmaSimBoot object.
#'
#'@param x An object of the GarmaSimBoot class, as provided
#'by \code{\link{GarmaSimBoot}}.
#'


setMethod(f="print",
          signature = "GarmaSimBoot",
          definition = function(x) {
            cat(" -------------------------------------------------------\n",
                "A",paste(x@sim@spec@family,"-Garma(",
                          x@sim@order[1],",",x@sim@order[2],")",
                          sep=""),"simulation bootstrap object: \n\n",

                if(x@sim@order[1] != 0) paste(
                  if(x@sim@order[1] == 1) "phi" else
                    paste("phi",1:x@sim@order[1],
                          " = ",sep=""), x@sim@spec@phi,"\n"),

                if(x@sim@order[2] != 0) paste(
                  if(x@sim@order[2] == 1) "theta" else
                    paste("theta",1:x@sim@order[2],
                          " = ",sep=""), x@sim@spec@theta,"\n"),

                if(TRUE%in%(colnames(x@sim@spec@X) != "null.vector"))
                  paste(paste("beta.",colnames(x@sim@spec@X), " = ",
                              sep=""), x@sim@spec@beta.x,"\n"),

                if(TRUE%in%(x@sim@spec@mu0 != 0)) paste(paste("mu0[", #altersr
                                                              1:length(x@sim@spec@mu0),"] = ",sep=""),
                                                        x@sim@spec@mu0,"\n"),

                if((x@sim@spec@family == "GA")&&(TRUE%in%(x@sim@spec@sigma2 != 0)))
                  paste("sigma2 = ", paste(x@sim@spec@sigma2),"\n"),

                "Number of Monte Carlo Simulations ('nmonte') = ",x@sim@nmonte,"\n",
                "Time Series Length ('nsteps') = ",x@sim@nsteps,"\n",
                "Burn-in ('burnin') = ",x@sim@burnin,"\n",
                "parallel = ",x@allow.parallel,"\n",
                "errorhandling = ",x@errorhandling,"\n",
                "n.try = ",x@n.try,"\n",
                "block length = ", x@l,"\n",
                "R = ",x@R,"\n",

                "-------------------------------------------------------\n",
                "\n",
                "Estimated parameters: \n\n")

            x@print.out

          }
)




#'@describeIn GarmaSimBoot
#'
#' Summary method for a GarmaSimBoot object.
#'
#'@param object An object of the GarmaSimBoot class, as provided
#'by \code{\link{GarmaSimBoot}}.
#'

setMethod(f="summary",
          signature = "GarmaSimBoot",
          definition = function(object,...) {

            cat(" -------------------------------------------------------\n",
                "A",paste(object@sim@spec@family,"-Garma(",
                          object@sim@order[1],",",object@sim@order[2],") ",
                          sep=""),"bootstrap object: \n\n",
                "Number of Monte Carlo Simulations ('nmonte') = ",object@sim@nmonte,"\n",
                "Time Series Length ('nsteps') = ",object@sim@nsteps,"\n",
                "Burn-in ('burnin') = ",object@sim@burnin,"\n",
                "parallel = ",object@allow.parallel,"\n",
                "errorhandling = ",object@errorhandling,"\n",
                "n.try = ",object@n.try,"\n",
                "block length = ", object@l,"\n",
                "R = ",object@R,"\n",

                "-------------------------------------------------------\n",
                "\n",
                "Bootstrap Statistic Summary: \n\n")

            if(object@sim@nmonte > 1)
            {
              db <- ddply(object@plot.out,.variables=c('.id','parameter','variable'), summarise, mean = mean(value))
              names(db)[1] <- 'length'
              db2 <- lapply(levels(db$variable), function(j) subset(db,variable==j))
              names(db2) <- levels(db$variable)

              return(db2)
            }
          })







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
          signature = "GarmaSimBoot",
          definition = function(x,scales = "free",...) {
            if(!scales%in%c('fixed', 'free_x', 'free_y', 'free'))
              stop("invalid value of scales, accepted values:
                   fixed', 'free_x', 'free_y', 'free' ")

            db <- x@plot.out[['db']]
            db <- db[db$variable == "original",]
            db2 <- x@plot.out[['db2']]

            if((length(x@l) > 1)&&(x@sim@nmonte > 1))
            {
              names(db2)[1] <- "true.value"

              if((length(x@sim@spec@beta.x) > 1)||
                 (x@sim@order[1] >1)||(x@sim@order[2] >1)) {
                ggplot(db, aes(value, fill = parameter)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  facet_grid(length~label,scales=scales) +
                  ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                   x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                " Average Bootstrap Estimates from ",x@sim@nmonte,
                                "Monte Carlo Simulation")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5))

              } else {
                ggplot(db, aes(value, fill = parameter)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  facet_grid(length~.,scales=scales) +
                  ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                   x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                " Average Bootstrap Estimates from ",x@sim@nmonte,
                                "Monte Carlo Simulations")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5))

              }

            } else {
              if(x@sim@nmonte > 1) {

                if((length(x@sim@spec@beta.x) > 1)||
                   (x@sim@order[1] >1)||(x@sim@order[2] >1)) {
                  ggplot(db, aes(value, fill = parameter)) +
                    geom_density(alpha = 0.2) +
                    geom_vline(data=db2,aes(xintercept=value,
                                            colour=true.value),size=1,linetype="dashed")+
                    facet_grid(label,scales=scales) +
                    ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                     x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                  " Average Bootstrap Estimates from ",x@sim@nmonte,
                                  "Monte Carlo Simulations")) +
                    theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                  face="bold.italic", colour="gray26",hjust=0.5))

                } else {
                  ggplot(db, aes(value, fill = parameter)) +
                    geom_density(alpha = 0.2) +
                    geom_vline(data=db2,aes(xintercept=value,
                                            colour=true.value),size=1,linetype="dashed")+
                    ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                     x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                  " Average Bootstrap Estimates from ",x@sim@nmonte,
                                  "Monte Carlo Simulations")) +
                    theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                  face="bold.italic", colour="gray26",hjust=0.5))

                }




              } else  {
                names(db2)[1] <- "true.value"
                db <- x@value[[1]][[1]]$t
                colnames(db) <- names(x@value[[1]][[1]]$t0)
                db <- reshape2::melt(db)[-1]
                names(db)[1] <- "parameter"

                ggplot(db, aes(value, fill = parameter)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                   x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                " Bootstrap Distribution from ",x@sim@nmonte,
                                "Monte Carlo Simulation")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5)) }


            }

          }
)

