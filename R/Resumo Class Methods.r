#### Spec
setClass("GarmaSpec",
	slots = list( beta.x = "numeric",
				 phi = "numeric", theta = "numeric",
				 X = "matrix"),
	prototype = list(beta.x = 1L,
				phi = 0L, theta = 0L),
	validity = function(object) {
	
	return(TRUE)				}
	
								)
								
								
							
setClass("PoissonSpec",
	slots = list(family = "character",  alpha = "numeric",
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
								

								
setClass("GammaSpec",
	slots = list(family = "character", sigma2 =  "numeric",
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
								
								
setGeneric("GarmaSpec", function(family, ...) {
	standardGeneric("GarmaSpec")
	})

	
	
setMethod(f="GarmaSpec",

	definition = function(family, ...) {
		
		switch(family,"PO" = {
			new("PoissonSpec",family = family, ...) 

							 }
				, "GA" = {
			new("GammaSpec",family = family, ...) 				
				
				         }						 
				)				
										} 
	
		)
		


#### Simulation		
		
setClass("GarmaSim",
	slots = list(spec = "GarmaSpec", nmonte = "numeric", 
		nsteps = "numeric", burnin = "numeric",
		allow.parallel = 'logical', seed = "numeric",
		 value = "list", order = "numeric"),
	prototype = list(nmonte = 1000, nsteps = 100,
		burnin = 1000, allow.parallel = TRUE,
		seed = 123),
	validity = function(object) {	
	if(!validObject(object@spec))
		return("spec must be a valid GarmaSpec object")	
		
	if(object@nmonte < 1)
		return("The number of Monte Carlo simulations
		must be a positive integer")
		
	if(object@nsteps < 1)
		return("The number of steps (length of the 
		simulated series) must be a positive integer")
		
	if(object@burnin < 0)
		return("The burnin must be zero or a positive integer")		
		
	return(TRUE)
								}						
								)	


setGeneric("GarmaSim", function(spec, ...) {
	standardGeneric("GarmaSim")
	})

	
	
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
	db[seq.int(from = nrow(db)-obj@nsteps + 1,
		to = nrow(db)),]
			
										  }
	slot(obj,"value") <- out
	
	return(obj)
									} 	
		)	
	

		
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
		{db <- sapply(x@value, function(j) j$yt)
		colnames(db) <- paste('yt',seq_len(x@nmonte),sep="_")
		rownames(db) <- seq_len(x@nsteps)
		db
		} else
			x@value[[1]]$yt
	
							}
		)
	
	
setMethod(f="plot",
	signature = "GarmaSim",
	definition = function(x,confInt=0.95,...) {
		
	if(x@nmonte > 1)
		{db <- sapply(x@value, function(j) j$yt)
		colnames(db) <- paste('yt',seq_len(x@nmonte),sep="_")
		rownames(db) <- seq_len(x@nsteps)
		db <- cbind(data.frame(index_t = seq_len(x@nsteps), mean_yt = 
		apply(db,1,mean)),as.data.frame(t(apply(db,1,quantile,
		probs=c((1-confInt)/2,0.5,confInt+(1-confInt)/2),type=1))))
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
		{db <- sapply(object@value, function(j) c(min(j$yt),mean(j$yt),max(j$yt)))
		db1 <- as.data.frame(t(apply(db,1,quantile,
		probs=c((1-0.95)/2,0.5,0.95+(1-0.95)/2),type=1)))
		db1 <- round(db1,4)
		db <- round(apply(db,1,mean),4)
		rownames(db1) <- NULL
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
			
			
### FIT

setClass("GarmaFit",
	slots = list(garma = 'GarmaSim', 
		allow.parallel = 'logical', seed = 'numeric',
		errorhandling = 'character' , n.try = 'numeric', 
		control = 'list',print.out = 'matrix', 
		plot.out = 'list', value = 'list'),
	prototype = list(allow.parallel = TRUE,
		seed = 123, errorhandling = 'try',
		n.try = 5L, control = list(iter.max=1000)),
	validity = function(object) {	
	if(!object@errorhandling %in% c("try","pass"))
		return("unrecognised value of 'errorhandling',
		accepted values are: 'try' and 'pass'")	
	
	if((object@errorhandling == "try")&&(object@n.try <= 0))
		return("'n.try' must be a positive integer")		
	
	if(object@allow.parallel)
		if (foreach::getDoParRegistered()==FALSE)
			return("parallel backend must be registered")
	return(TRUE)
								}						
								
		)		

		
setGeneric("GarmaFit", function(garma, ...) {
	standardGeneric("GarmaFit")
	})

	

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
					obj@garma@spec@X )),
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
					obj@garma@spec@X )),
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


		
setMethod(f="plot",
	signature = "GarmaFit",
	definition = function(x,scales = "free",...) {
	if(!scales%in%c('fixed', 'free_x', 'free_y', 'free'))
		stop("invalid value of scales, accepted values: 
		fixed', 'free_x', 'free_y', 'free' ")
	
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
		"parameters",summarise,Min. = min(value),
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
		
		
		
###BOOT

setClass("GarmaSimBoot",
	slots = list(sim = 'GarmaSim', l = 'numeric',
		R = 'numeric', allow.parallel = 'logical',
		seed = 'numeric', errorhandling = 'character',
		n.try = 'numeric', boot.function = 'function',		
		control = 'list',print.out = 'data.frame', 
		plot.out = 'list', value = 'list'),
	prototype = list( R = 100L, l = integer(0),
		allow.parallel = TRUE,seed = 123, 
		errorhandling = 'try', n.try = 5L,
		boot.function = function(x) 0,
		control = list(iter.max=1000)),
	validity = function(object) {	
	if(identical(object@l,integer(0)))
		return("'l' must be a positive integer smaller
		than nsteps.")
		
	if((object@l <= 0)||(object@l > object@sim@nsteps))
		return("'l' must be a positive integer smaller
		than nsteps.")

	if(object@R <=0)
		return("'R' must be a positive integer.")

	if(!object@errorhandling %in% c("try","pass"))
		return("unrecognised value of 'errorhandling',
		accepted values are: 'try' and 'pass'.")	
	
	if((object@errorhandling == "try")&&(object@n.try <= 0))
		return("'n.try' must be a positive integer.")		
	
	if(object@allow.parallel)
		if (foreach::getDoParRegistered()==FALSE)
			return("parallel backend must be registered.")
	return(TRUE)
								}						
								
		)		

							
setGeneric("GarmaSimBoot", function(sim, ...) {
	standardGeneric("GarmaSimBoot")
	})



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
							obj@sim@spec@X )),
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
						obj@sim@spec@X )), 
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
	
	clean <- llply(obj@value,function(j) {
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
		melt(i,id.vars="parameter")
		})})
		
	slot(obj,"plot.out") <- plot.out
	

	return(obj)							   
							  
									}
	)



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
			

		