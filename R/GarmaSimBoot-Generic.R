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
#'


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
#'
#'@examples ##Some Specs/Simulations for different outputs:
#'Sim1 <- GarmaSim(
#'GarmaSpec("GA",
#'phi = 0.5,
#'beta.x = 1,
#'sigma2 =2,
#'X = as.matrix(
#'  data.frame(
#'    x1 = rep(10,101)))), 
#'nmonte = 1, burnin = 0)
#'
#'
#'Sim2 <- GarmaSim(
#'GarmaSpec("PO",
#'phi = 0.2,
#'theta = c(0.1, 0.3, 0.5),
#'mu0 = 1:3), 
#'nmonte = 10, burnin = 0)
#'
#'
#'Sim3 <- GarmaSim(
#'GarmaSpec("PO",
#'beta.x = c(1,1),
#'X = as.matrix(
#'  data.frame(
#'  intercept = rep(1,1100),
#'    x1 = c(rep(7,100),
#'        rep(2,1000))))), 
#'nmonte = 10)
#'
#'
#'Sim4 <- GarmaSim(
#'GarmaSpec("GA",
#'phi = 0.5,
#'beta.x = 1,
#'sigma2 =2,
#'X = as.matrix(
#'  data.frame(
#'    x1 = rep(10,101)))), 
#'nmonte = 10, burnin = 0)
#'
#'
#' ##Example of the GarmaSimBoot methods:
#'
#'# one monte carlo, one block length:
#'ex1 <- GarmaSimBoot(
#'Sim1, l = 20)
#'print(ex1)
#'plot(ex1)
#'summary(ex1)
#'
#'
#'# one monte carlo, multiple block lengths +
#'# a user defined function to apply in the 
#'# bootstrap resamples (0.1 and 0.9 quantiles):
#'ex2 <- GarmaSimBoot(
#'Sim1, 
#'l = c(4,7,10),
#'boot.function = 
#'  function(x) quantile(x,
#'     probs = c(0.1,0.9))
#')
#'print(ex2)
#'plot(ex2)
#'summary(ex2)
#'
#'          
#'# 10 monte carlo sim, multiple block lengths +
#'# max(order(garma)) > 1:
#'ex3 <- GarmaSimBoot(
#'Sim2, 
#'l = c(10,15))
#'print(ex3)
#'plot(ex3)
#'plot(ex3, variable = "Median.")
#'plot(ex3, type = "original-bias")
#'summary(ex3)
#'
#'
#'# 10 monte carlo sim, multiple block lengths +
#'# order = c(0,0); no AR/MA term:
#'ex4 <- GarmaSimBoot(
#'Sim3, 
#'l = c(15,20,30))
#'print(ex4)
#'plot(ex4)
#'plot(ex4, variable = "1st Qu.")
#'plot(ex4, type = "original-bias")
#'summary(ex4)
#'
#'
#'# 10 monte carlo sim, single block lengths +
#'# ARMA+Beta.x
#'ex5 <- GarmaSimBoot(
#'Sim4, 
#'l = c(20))
#'print(ex5)
#'plot(ex5)
#'plot(ex5, type = "original-bias")
#'summary(ex5)
#'
#'
#'# 10 monte carlo sim, two block lengths +
#'# ARMA+Beta.x
#'ex6 <- GarmaSimBoot(
#'Sim4, 
#'l = c(5,30))
#'print(ex6)
#'plot(ex6)
#'plot(ex6, type = "original-bias")
#'summary(ex6)
#'





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

            clean <- if(obj@sim@nmonte > 1)
              plyr::llply(obj@value,function(j) {
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
            })else
              plyr::llply(obj@value,function(j) {
                boot.function <- obj@boot.function
                plyr::llply(j,.parallel = obj@allow.parallel,
                            .fun = function(x) {
                              with(x,{
                                colnames(t) <- names(t0)
                                t
                              })
                            })
              }
              )

            print.out <- {

              out <- if(obj@sim@nmonte > 1) {
                out <- plyr::ldply(clean,function(j) {
                  plyr::ldply(j,
                              function(i){
                                i$parameter <- rownames(i)
                                i })})
                names(out)[1] <- "length"
                out <- out[c(1,ncol(out),seq.int(2,ncol(out)-1))]
                out
                                            } else {
                  out <- plyr::ldply(clean, reshape2::melt)[,-c(2,5)]
                  colnames(out) <- c("length","parameter","original")
                  out
                      }

              if("beta.null.vector"%in%out$parameter) {
                ind <-  "beta.null.vector" != out$parameter
                out <- out[ind,]
                                                      }

              if(all(obj@sim@order == c(0,0)))
                out$parameter <- paste("beta.",out$parameter,sep="")

              as.data.frame(out)
            }
            slot(obj,"print.out") <- print.out


            plot.out  <- {

              create.labels <- function(data,labels =
                                          c("beta","phi","theta","sigma")) {
                data$label <- 0
                for (i in labels) {
                  data$label[grep(i,data$parameter)] <- i
                }
                return(data)
              }

              create.parameters <- function(x) {
                par1 <- if(identical(grep("beta", x@print.out$parameter),
                                     integer(0))) NULL else x@sim@spec@beta.x
                par2 <- if(identical(grep("phi", x@print.out$parameter),
                                     integer(0))) NULL else x@sim@spec@phi
                par3 <- if(identical(grep("theta", x@print.out$parameter),
                                     integer(0))) NULL else x@sim@spec@theta
                par4 <- if(identical(grep("sigma", x@print.out$parameter),
                                     integer(0))) NULL else sqrt(x@sim@spec@sigma2)
                pars <- c(par1,par2,par3,par4)
                return(pars)
              }

              db <- if(obj@sim@nmonte > 1) {
                plyr::ldply(clean,
                            function(j) {
                              plyr::ldply(j,
                                          function(i){
                                            i$parameter <- rownames(i)
                                            reshape2::melt(i,id.vars="parameter")
                                          })})} else {
                        out <- plyr::ldply(clean, reshape2::melt)[,-c(2,5)]
                        colnames(out) <- c("length","parameter","value")
                        out$variable <- factor("original")
                        out

                                          }
              if("beta.null.vector"%in%db$parameter) {
                ind <-  "beta.null.vector" != db$parameter
                db <- db[ind,]
              }

              if(all(obj@sim@order == c(0,0)))
                db$parameter <- paste("beta.",db$parameter,sep="")
              colnames(db)[1] <-"length"
              db <- create.labels(db)
              db2 <- data.frame(parameter = unique(db$parameter),
                               variable = factor("true.value"),
                               value = create.parameters(obj))
              db2 <- create.labels(db2)



              list(db = db, db2=db2) #db2



            }

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
              db <- plyr::ddply(
                object@plot.out$db,
                .variables = c('length','parameter','variable'),
                'summarise',
                mean = mean(value))


              colnames(db)[1] <- 'length'
              db2 <- lapply(levels(db$variable), function(j) subset(db,variable==j))
              names(db2) <- levels(db$variable)
			  
			  db <- object@plot.out$db2[c('parameter','value')]
              names(db)[2] <- 'true.value'
			
              return(list("True values"= db,
				Estimates = db2))
            } else {
              db2 <- object@plot.out$db2[c('parameter','value')]
              names(db2)[2] <- 'true.value'
			  db <- plyr::ddply(
			  object@plot.out$db, 
			  .variables = c("length","parameter","variable"),
			  'summarise',
			  mean=mean(value))
			  colnames(db)[4] <- "Mean Est."	
              db <- db[order(db$length, decreasing=T),][order(db$variable),]
			  
			  out <- list("True values" = db2,
				Estimates = db)
              return(out)
            }

          })

#'@describeIn GarmaSimBoot
#'
#' Plot method for a GarmaSimBoot object.
#'
#'@param x An object of the GarmaSimBoot class, as provided
#'by \code{\link{GarmaSimBoot}}.
#'
#'@param scales A character specifying if the scales should be free or not.
#'The default value is \code{"free"} and accepted values are:
#' \code{'fixed', 'free_x', 'free_y' and 'free'}.
#' 
#' @param variable A string with the name of the statistic to plot.
#' The default value is the Mean of the replicates. Also, note
#' that if \code{x@sim@nmonte = 1} there is no statistic being 
#' computed and only the default value of
#' 'original' is accepted. Any other option for \code{variable}
#' is overwritten. 
#' 
#' @param type A string with the value 'default' 
#' or 'original-bias' designating the type of the plot. The 
#' 'default' plot is the original statistic, while the 
#' 'original-bias' contrast the original with the bias-corrected
#'  estimates.

setMethod(f="plot",
          signature = "GarmaSimBoot",
          definition = function(x,scales = "free",
                                variable= "Mean.",
                                type = "default",...) {
            if(!scales%in%c('fixed', 'free_x', 'free_y', 'free'))
              stop("invalid value of scales, accepted values:
                   fixed', 'free_x', 'free_y', 'free' ")
            if(x@sim@nmonte == 1)
            {
              variable <- "original"
              type = "default"
            }
            

            db <- x@plot.out[['db']]

            if(!(variable%in%levels(db$variable)))
              stop(paste("Accepted values of 'variable' are:",levels(db$variable)))
            
            if(!type%in%c("default","original-bias"))
              stop("Accepted values of 'type' are:
                         'default' and 'original-bias'")
            if(type == "original-bias")
              db <- db[(db$variable=='original')|
                         (db$variable=='bias corrected'),] else 
              db <- db[db$variable == variable,]         
                           
                         
            
            ind <- !(is.na(db$value)|is.nan(db$value)|is.infinite(db$value))
            db <- db[ind,]
            db2 <- x@plot.out[['db2']]
            db2 <-db2[db2$parameter%in%unique(db$parameter),]
            
            
            if(type == 'default')
            { names(db2)[1] <- "true.value"
              if((length(x@l) > 1)&&(x@sim@nmonte > 1))
            {
              
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
                                "Monte Carlo Simulations")) +
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
                    facet_grid(label~.,scales=scales) +
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
                    facet_grid(length~.,scales=scales)+
                    ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                     x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                  " Average Bootstrap Estimates from ",x@sim@nmonte,
                                  "Monte Carlo Simulations")) +
                    theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                  face="bold.italic", colour="gray26",hjust=0.5))

                }




              } else  {

                ggplot(db, aes(value, fill = parameter)) +
                  geom_density(alpha = 0.2) +
                  geom_vline(data=db2,aes(xintercept=value,
                                          colour=true.value),size=1,linetype="dashed")+
                  facet_grid(length~.,scales=scales)+
                  ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                   x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                                " Bootstrap Distribution from ",x@sim@nmonte,
                                "Monte Carlo Simulation")) +
                  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                                face="bold.italic", colour="gray26",hjust=0.5)) }


            }

            } else {
              db3 <- plyr::ddply(db,
                                .variables = c("length","parameter","variable"),
                                'summarise',
                                mean=mean(value))
              db3 <- db3[db3$variable=='original'|db3$variable=='bias corrected',]
              colnames(db3)[4] <- "value"
              db4 <- do.call("rbind", replicate(length(unique(db$length)), db2[,-4], simplify = FALSE)) 
              db4$length <- sort(rep(unique(db$length),nrow(db2)))

              db2 <- rbind(db3,db4)
              db2$variable <- paste("Mean - ",db2$variable, sep = "")
                
              ggplot(db, aes(value, fill = variable)) +
                geom_density(alpha = 0.2) +
                geom_vline(data=db2,aes(xintercept=value,
                                        colour=variable),size=1,linetype="dashed")+
                facet_grid(length~parameter,scales=scales) +
                ggtitle(paste("Simulated", paste(x@sim@spec@family,"-Garma(",
                                                 x@sim@order[1],",",x@sim@order[2],"): ",	sep=""),
                              " Average Bootstrap Estimates from ",x@sim@nmonte,
                              "Monte Carlo Simulations")) +
                theme(plot.title=element_text(size=rel(1.2), lineheight=.9,
                                              face="bold.italic", colour="gray26",hjust=0.5))
                    }
            
          }
)




