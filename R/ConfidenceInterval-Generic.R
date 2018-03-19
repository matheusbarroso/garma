#' @include ConfidenceInterval.R
NULL


#' ConfidenceInterval  class generic function.
#'
#' A generic function working as a class constructor for 
#' the ConfidenceInterva class.
#'
#' This generic function is designed to be a constructor for the
#' ConfidenceInterva class.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@param garma.boot  An object of the \bold{GarmaSimBoot} 
#'class, as provided by \code{\link{GarmaSimBoot}}.
#'
#'@param \dots Further arguments that specify the 
#'ConfidenceInterval object, check the \bold{slots}
#' arguments.
#'


setGeneric("ConfidenceInterval", function(garma.boot, ...) {
  standardGeneric("ConfidenceInterval")
})


#'@describeIn ConfidenceInterval
#'
#' ConfidenceInterval Class constructor.
#'
#' A method to create an object of the ConfidenceInterval class.
#'
#' This method creates an object of the ConfidenceInterval class.
#'
#'@param garma.boot  An object of the \bold{GarmaSimBoot} 
#'class, as provided by \code{\link{GarmaSimBoot}}.
#'
#'@param \dots Further arguments that specify the Confidence
#'Interval object,check the \bold{slots} arguments.
#'
#'@examples ##Some Specs/Simulations for different boot outputs:
#'library(garma)
#'no_cores <- if(detectCores()==1) 1 else detectCores() -1
#'registerDoParallel(no_cores)
#'
#'spec1 <- GarmaSim(
#'GarmaSpec("GA",
#'          phi = c(0.5),
#'          beta.x = 1,
#'          sigma2 = 2,
#'          X = as.matrix(
#'            data.frame(
#'              x1 = rep(10,101)))),
#'nmonte = 1, burnin = 0)
#'
#'boot1 <- GarmaSimBoot(spec1,l = 20 )
#'
#'spec2 <- GarmaSim(
#'GarmaSpec("GA",
#'          phi = c(0.5),
#'          beta.x = 1,
#'          sigma2 = 2,
#'          X = as.matrix(
#'            data.frame(
#'              x1 = rep(10,101)))),
#'nmonte = 10, burnin = 0)
#'
#'boot2 <- GarmaSimBoot(spec2,l = c(5,15,20) )
#'
#' ##Example of the GarmaSimBoot methods:
#'
#'ci <- ConfidenceInterval(boot1)
#'print(ci)
#'summary(ci)
#'coverage (ci)
#'plot(ci,type.ci='basic')
#'
#'ci <- ConfidenceInterval(boot2)
#'print(ci)
#'summary(ci)
#'coverage (ci)
#'plot(ci,type.ci='perc')
#'plot(ci)  #norm ci


setMethod(f="ConfidenceInterval",
          definition = function(garma.boot, ...) {
            obj <- new("ConfidenceInterval",
                       garma.boot = garma.boot, ...)
            library(plyr)
            
            conf <- obj@conf
            
            out <- 	{
              
              intervals <- lapply(obj@garma.boot@value, function(j){
                conf <- conf
                llply(j, 
                      .parallel = obj@allow.parallel,
                      #.paropts = list(.export = c('conf')),
                      .fun = function(i) {
                        t0 <- i$t0
                        tn <- i$t
                        npar <- ncol(tn)
                        ldply(seq_len(npar), 
                              function (k) {
                                
                                bias.c <- boot:::norm.ci(
                                  t0 = t0[k],
                                  t = tn[,k], 
                                  conf = conf)
                                
                                basic <- boot:::basic.ci(
                                  t0 = t0[k],
                                  t = tn[,k], 
                                  conf = conf)
                                
                                perc <- boot:::perc.ci(
                                  t = tn[,k],
                                  conf = conf)
                                
                                data.frame(bias.c.LB = bias.c[2],
                                           bias.c.UB = bias.c[3],
                                           basic.LB = basic[4],
                                           basic.UB = basic[5],
                                           perc.LB = perc[4],
                                           perc.UB = perc[5])
                                })
                        
                                        })
                })
              
              fit <- GarmaFit(obj@garma.boot@sim)@value
              Z.lb <- qnorm((1-conf)/2)
              Z.ub <- qnorm((1-conf)/2+conf)
              norm <- llply(fit, 
                            .parallel = obj@allow.parallel, 
                            #.paropts = list(.export = c('Z.lb','Z.ub')),
                            .fun = function(j)
                              data.frame( 
                                norm.LB = j$coef + sqrt(diag(j$vcov))*Z.lb,
                                norm.UB = j$coef + sqrt(diag(j$vcov))*Z.ub))
              
              cis <- llply(seq_len(length(obj@garma.boot@l)), 
                           function(i) {
                             norm <- norm
                             #conf.interv <- conf.interv
                             ldply(seq_len(obj@garma.boot@sim@nmonte),
                                   .parallel = obj@allow.parallel,
                                  # .paropts = list(.export = c('norm','conf.interv')),
                                   .fun = function(j) {
                                     
                                     conf.interv <- cbind(
                                       norm[[j]],
                                       intervals[[i]][[j]])
                                     
                                     conf.interv$conf <- conf
                                     conf.interv$parameter <- rownames(norm[[j]])
                                     conf.interv$nmonte <- j
                                     conf.interv
                                     })
                           })
              
              names(cis) <- paste("l=",obj@garma.boot@l,sep="")
              
              cis
              
              }
            
            slot(obj,"value") <- out # list
            

            plot.out  <- {
            
              db <- ldply(obj@value,
                          reshape2::melt, 
                          id.vars=c("parameter","nmonte"))  
              names(db)[1] <- "length"
              db
            }
            
            slot(obj,"plot.out") <- plot.out #data.frame
            
            summary.out  <- {
              
              ddply(obj@plot.out,
                    .variables = c('length','parameter','variable'), 
                    summarise, 
                    mean = mean(value)
                    )
              
            }
            
            slot(obj,"summary.out") <- summary.out #data.frame
            
            
            
            
            coverage.out  <- {
              
              pars <- unique(ldply(obj@value,
                                   reshape2::melt, 
                                   id.vars=c("parameter","nmonte")
                                   )$parameter)
              ## phi ##
              phi <- 
                if(obj@garma.boot@sim@order[1] == 0) 0L else 
                  if(obj@garma.boot@sim@order[1] == 1)
                    data.frame(phi = obj@garma.boot@sim@spec@phi) else
                      {db <- data.frame(t(obj@garma.boot@sim@spec@phi))
                      names(db) <- paste("phi",seq_len(obj@garma.boot@sim@order[1]),sep = "")
                      }
              
              ## theta ##
              theta <- 
                if(obj@garma.boot@sim@order[2] == 0) 0L else
                  if(obj@garma.boot@sim@order[2] == 1)
                    data.frame(theta = obj@garma.boot@sim@spec@theta) else
                      {db <- data.frame(t(obj@garma.boot@sim@spec@theta))
                      names(db) <- paste("theta",seq_len(obj@garma.boot@sim@order[2]),sep = "")
                      }
              
              
              ## sigma2 ##
              sigma <- 
                if(obj@garma.boot@sim@spec@family == "GA")
                  data.frame(sigma = sqrt(obj@garma.boot@sim@spec@sigma2)) else 0L
              
              
              ## beta ##		
              beta.x <-
                {db <- data.frame(t(obj@garma.boot@sim@spec@beta.x))
                names(db) <- paste("beta.",colnames(obj@garma.boot@sim@spec@X), sep = "")
                db
                }
              
              true.pars <- cbind(phi,theta,sigma,beta.x)
              
              
              coverage <- 
                ldply(obj@value, .fun = function(i) {
                  ldply(pars, .fun = function(j){
                    
                    true.value <- true.pars[[j]]
                    
                    db <- subset(i, parameter == j)
                    
                    norm <- (db$norm.LB <= true.value)&(db$norm.UB >= true.value) 
                    bias.c <- (db$bias.c.LB <= true.value)&(db$bias.c.UB >= true.value)
                    basic <- (db$basic.LB <= true.value)&(db$basic.UB >= true.value)
                    perc <-  (db$perc.LB <= true.value)&(db$perc.UB >= true.value)
                    
                    coverage <- sapply(list(norm,bias.c,basic,perc),sum)/sapply(list(norm,bias.c,basic,perc),length)
                    names(coverage) <- c('norm','bias.c','basic','perc')
                    coverage <- cbind(as.data.frame(t(coverage)),parameter=j)
                    names(coverage) <- c('norm','bias.c','basic','perc','parameter')
                    
                    coverage
                    })
                  })
              
              
              names(coverage)[1] <- "length"
              coverage[,c(1,6,2:5)]
              
            }
            
            slot(obj,"coverage.out") <- coverage.out #data.frame
            
            
            return(obj)
            
          }
)


#'@describeIn ConfidenceInterval
#'
#' Print method for a ConfidenceInterval object.
#'
#'@param x An object of the ConfidenceInterval class, as provided
#'by \code{\link{ConfidenceInterval}}.
#'

setMethod(f="print",
          signature = "ConfidenceInterval",
          definition = function(x) {
            cat(" -------------------------------------------------------\n",
                "A",paste(x@garma.boot@sim@spec@family,"-Garma(",
                          x@garma.boot@sim@order[1],",",x@garma.boot@sim@order[2],")",
                          sep=""),"simulation bootstrap object: \n\n",
                
                if(x@garma.boot@sim@order[1] != 0) paste(
                  if(x@garma.boot@sim@order[1] == 1) "phi" else
                    paste("phi",1:x@garma.boot@sim@order[1],
                          " = ",sep=""), x@garma.boot@sim@spec@phi,"\n"),
                
                if(x@garma.boot@sim@order[2] != 0) paste(
                  if(x@garma.boot@sim@order[2] == 1) "theta" else
                    paste("theta",1:x@garma.boot@sim@order[2],
                          " = ",sep=""), x@garma.boot@sim@spec@theta,"\n"),
                
                if(TRUE%in%(colnames(x@garma.boot@sim@spec@X) != "null.vector"))
                  paste(paste("beta.",colnames(x@garma.boot@sim@spec@X), " = ",
                              sep=""), x@garma.boot@sim@spec@beta.x,"\n"),
                
                if(TRUE%in%(x@garma.boot@sim@spec@mu0 != 0)) paste(paste("mu0[", #altersr
                                                              1:length(x@garma.boot@sim@spec@mu0),"] = ",sep=""),
                                                        x@garma.boot@sim@spec@mu0,"\n"),
                
                if((x@garma.boot@sim@spec@family == "GA")&&(TRUE%in%(x@garma.boot@sim@spec@sigma2 != 0)))
                  paste("sigma2 = ", paste(x@garma.boot@sim@spec@sigma2),"\n"),
                
                "Number of Monte Carlo Simulations ('nmonte') = ",x@garma.boot@sim@nmonte,"\n",
                "Time Series Length ('nsteps') = ",x@garma.boot@sim@nsteps,"\n",
                "Burn-in ('burnin') = ",x@garma.boot@sim@burnin,"\n",
                "parallel = ",x@allow.parallel,"\n",
                "block length = ", x@garma.boot@l,"\n",
                "R = ",x@garma.boot@R,"\n",
                
                "-------------------------------------------------------\n",
                "\n",
                "Estimated parameters: \n\n")
            
            x@value
            
            
          }
)

#'@describeIn ConfidenceInterval
#'
#' Summary method for a ConfidenceInterval object.
#'
#'@param object An object of the ConfidenceInterval class, as provided
#'by \code{\link{ConfidenceInterval}}.
#'

setMethod(f="summary",
          signature = "ConfidenceInterval",
          definition = function(object,...) {
            
            cat(" -------------------------------------------------------\n",
                "A",paste(object@garma.boot@sim@spec@family,"-Garma(",
                          object@garma.boot@sim@order[1],",",object@garma.boot@sim@order[2],") ",
                          sep=""),"bootstrap object: \n\n",
                "Number of Monte Carlo Simulations ('nmonte') = ",object@garma.boot@sim@nmonte,"\n",
                "Time Series Length ('nsteps') = ",object@garma.boot@sim@nsteps,"\n",
                "Burn-in ('burnin') = ",object@garma.boot@sim@burnin,"\n",
                "parallel = ",object@garma.boot@allow.parallel,"\n",
                "block length = ", object@garma.boot@l,"\n",
                "R = ",object@garma.boot@R,"\n",
                
                "-------------------------------------------------------\n",
                "\n",
                "Bootstrap Statistic Summary: \n\n")
            
            object@summary.out
          })

         

#'@describeIn ConfidenceInterval
#'
#' Plot method for a ConfidenceInterval object.
#'
#'@param x An object of the ConfidenceInterval class, as provided
#'by \code{\link{ConfidenceInterval}}.
#'
#'@param scales A character specifying if the scales should be free or not.
#'The default value is \code{"free"} and accepted values are:
#' \code{'fixed', 'free_x', 'free_y' and 'free'}.
#' 
#' 
#' @param type.ci A string designating the confidence interval to plot, only 
#' aplicable when \code{nmonte > 1}. The true value is 'norm' and the other accepted
#' values are: 'bias.c', 'basic' and 'perc'. Refering to the Normal, 
#' Bias Corrected, Basic and Percentile ci's.

setMethod(f="plot",
          signature = "ConfidenceInterval",
          definition = function(x,scales = "free",
                                type.ci = "norm") {
            
            if(x@garma.boot@sim@nmonte == 1) {
              true.pars <- x@garma.boot@plot.out$db2
              
              temp <- rdply(4,true.pars)
              temp[,3] <- NULL
              temp$variable <- 'norm'
              temp$variable[temp$.n == 2] <- 'bias.c'
              temp$variable[temp$.n == 3] <- 'basic'
              temp$variable[temp$.n == 4] <- 'perc'
              temp$.n <- NULL
              temp <- rdply(length(x@garma.boot@l),temp)
              temp$length <- "temp"
              for (i in seq_len(length(x@garma.boot@l))) {
                temp$length[temp$.n == i] <- paste("l=",x@garma.boot@l[i],sep="")
              }
              
              
              bounds <- ldply(unique(x@plot.out$length), .fun = function(k) {
                ldply(unique(x@plot.out$parameter), .fun =  function(i){
                  db <- subset(temp, (parameter == i)&(length == k))
                  db2 <- subset(x@plot.out, (parameter == i)&(length == k) )
                  
                  bounds <- ldply(db$variable, function(j){
                    lb <- paste(j,".LB", sep='')
                    ub <- paste(j,".UB", sep='')
                    db3 <- db2$value[match(c(lb,ub),db2$variable)]
                    names(db3) <- c("LB","UB")
                    db3
                    })
                  })
                }) 
              
              db <- cbind(temp,bounds)
              
              ggplot(db, aes(y=value,x=variable)) + 
                geom_point() + 
                geom_errorbar(aes(ymin=LB, ymax=UB, x=variable))+
                facet_grid(length ~ parameter,scales = scales) +
                ggtitle(" CI's Lower Bounds (LB) and Upper Bounds 
                        (LB) for parameters true values") +
                theme(plot.title = element_text(size = rel(1.2), 
                                                lineheight=.9,
                                                face="bold.italic",
                                                colour="gray26",
                                                hjust=0.5))
            
                
            } else {
              if(!type.ci %in% c('norm','bias.c','basic','perc'))
                stop("Invalid type.ci value, accepted values are:
                     'norm','bias.c','basic','perc'")
              
              db <- x@plot.out
              db <- subset(db, 
                           ({variable == paste(type.ci, ".LB", sep="")}|
                              {variable == paste(type.ci, ".UB", sep="")}))
              
              db2 <- x@garma.boot@plot.out$db2

              ggplot(db, aes(value, fill = variable)) +
                geom_density(alpha = 0.2) +
                geom_vline(data = db2,
                           aes(xintercept=value),size = 1,linetype = "dashed")+
                facet_grid(length ~ parameter,scales = scales) + 
                ggtitle(paste(type.ci,
                              " CI's Lower Bounds (LB) and Upper Bounds (LB) for 'nsims' =",
                              x@garma.boot@sim@nmonte)) +
                theme(plot.title = element_text(size = rel(1.2), 
                                                lineheight=.9,
                                                face="bold.italic",
                                                colour="gray26",
                                                hjust=0.5))
            }
            }
          )


#@describeIn ConfidenceInterval
#'
#' coverage generic function.
#'
#' A generic function working to compute ci's coverage rates
#'
#' This generic function is designed to extract the coverage
#' rates from a \code{\link{ConfidenceInterval}} object.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@param object  An object of the \bold{ConfidenceInterval} 
#'class, as provided by \code{\link{ConfidenceInterval}}.
#'




setGeneric("coverage", function(object) {
  standardGeneric("coverage")
})

#'@describeIn ConfidenceInterval
#'
#' coverage method for a ConfidenceInterval object.
#'
#'@param object An object of the ConfidenceInterval class, as provided
#'by \code{\link{ConfidenceInterval}}.
#'

setMethod(f="coverage",
          signature = "ConfidenceInterval",
          definition = function(object) {
            
            cat(" -------------------------------------------------------\n",
                "A",paste(object@garma.boot@sim@spec@family,"-Garma(",
                          object@garma.boot@sim@order[1],",",object@garma.boot@sim@order[2],") ",
                          sep="")," MBB ci's coverage object: \n\n",
                "Number of Monte Carlo Simulations ('nmonte') = ",object@garma.boot@sim@nmonte,"\n",
                "Time Series Length ('nsteps') = ",object@garma.boot@sim@nsteps,"\n",
                "parallel = ",object@garma.boot@allow.parallel,"\n",
                "block length = ", object@garma.boot@l,"\n",
                "R = ",object@garma.boot@R,"\n",
                
                "-------------------------------------------------------\n",
                "\n",
                "Bootstrap Confidence Interval Coverage Rates: \n\n")
            
            object@coverage.out
          })





