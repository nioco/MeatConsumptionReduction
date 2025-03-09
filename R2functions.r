#Supplemental R Functions for R-square Change in SEMs with Latent Variables and Missing Data

rsquareCalc <- function(model, y, x, adj = FALSE, effN = FALSE, silent = FALSE){
	#model is a model fit by lavaan using e.g., the sem() or lavaan() function.
	#y is a character vector of length 1 specifying the name of the (single) structural outcome of interest.
	#x is a vector of one or more character strings specifying the name(s) of the target predictor(s) of interest, to be omitted from the reduced model.
	#adj: do you want to calculate adjusted rather than unadjusted R2 and R2 change? Defaults to FALSE.
	#effN: if TRUE, N in the adjusted R-square calculation is set to the lowest effective N in the structural regression. Defaults to FALSE.
	#silent: if TRUE, output does not automatically print (but is returned as invisible). Defaults to FALSE.
		#This argument is invoked when using rsquareCalc.Boot in order to turn off default printing while taking bootstrap resamples.
	
	require("lavaan")
	
	if(!is.character(y)|length(y) != 1) stop("y must be a character vector of length 1, specifying the name of the DV in the (manifest or latent variable) regression of interest!")
	
	#parameter estimates
	pe <- parameterEstimates(model, standardized = TRUE, rsquare = TRUE)
	
	#correlation matrix of all variables
	Rmat <- lavInspect(model, what = "cor.all")
	
	#regression coefficients
	Gamma <- pe[pe$lhs == y & pe$op == "~", ]
	
	#names of X variables NOT specified in x
	otherXnames <- Gamma[!(Gamma$rhs %in% x), "rhs"]

	
	#Grab correlation matrix of other Xs.
	Rxx <- Rmat[otherXnames, otherXnames, drop = FALSE]
	
	#Inverse X cor mat.
	RxxInv <- solve(Rxx)
	
	#vector of xy correlations.
	Rxy <- Rmat[otherXnames,y, drop = FALSE] #this way preserves the correct order of the other x names
	
	#compute new gammas as they would have been without the variables in x included in the model.
	#gamma = RxxInv%*Rxy
	GammaNew <- RxxInv%*%Rxy
	
	#R square of the submodel
	RsqReduced <- t(GammaNew)%*%Rxy
	
	#R square from Full model
	RsqFull <- pe[pe$lhs == y & pe$op == "r2", "est"]
	
	if(adj){
		#If adjusted R-square is requested
		#Retrieve number of observations used in the analysis.
		n <- lavInspect(model, what = "nobs")
		
		#Number of predictors in the full model.
		pFull <- nrow(pe[pe$lhs == y & pe$op == "~",])
		
		#Number of predictors contributing to increment in R-squared.
		pInc <- length(x)
		
		#Reducted model p = pFull - pInc.
		pRed <- pFull - pInc
		
		if(effN){
			
			#If effective N is requested, first check that the fmi is calculable.
			#To do this, the following code borrows from lavaan's internal code.
			
			###Code taken from parameterEstimates() function:
			PT <- parTable(model)
			EM.cov <- lavInspect(model, "sampstat.h1")$cov
            EM.mean <- lavInspect(model, "sampstat.h1")$mean
            this.options <- model@Options
            this.options$optim.method <- "none"
	        this.options$sample.cov.rescale <- FALSE
	        this.options$check.gradient <- FALSE
	        this.options$baseline <- FALSE
	        this.options$h1 <- FALSE
	        this.options$test <- FALSE
			fit.complete <- lavaan(model = PT, sample.cov = EM.cov, 
            sample.mean = EM.mean, sample.nobs = n, slotOptions = this.options)
            ###
            
            #Check that the complete model is identified:
            if(any(eigen(lavInspect(fit.complete, what = "vcov"))$values <0)){
           		
           		#If the model used to estimate the fmi is non-identified, everything is NA.
           		res <- rep(NA, 2)
           		names(res) <- c(paste("Rsquare Without ", paste0(x, collapse = " "), collapse = ""), "RsquareChange")
           		
            }else{        
	            
	            #Otherwise, the calculations proceed
	            
	            #peFMI = parameter estimates with fmi
	            peFMI <- parameterEstimates(model, standardized = TRUE, rsquare = TRUE, fmi = TRUE)
				
				#Flag regression relationships with proper dv:
				regressionflag <- peFMI$lhs == y & peFMI$op == "~"
				
				#Flag the residual variance as well
				residvarflag <- peFMI$lhs == y & peFMI$op == "~~" & peFMI$rhs == y
				
				#Subset the parameter estimates object to include all 
				#contributors to predicted and residual variance.
				peFMI_sub <- peFMI[regressionflag|residvarflag,]
				
				#Retrieve max fmi from structural model.
				fmi <- peFMI_sub[which.max(peFMI_sub$fmi), "fmi"]
				
				#Calculate the effective n.
				effN <- n*(1-fmi)
				
				#Overwrite the original n for use in subsequent calculations.
				n <- effN
				
				#Adjusted R-square calculations.
				multiplierFull <- (n-1)/(n - pFull - 1)
				multiplierRed <- (n-1)/(n - pRed - 1)
				RsqReduced <- 1 - multiplierRed*(1 - RsqReduced)
				RsqFull <- 1 - multiplierFull*(1 - RsqFull)
				
				#R square change is difference between overall R square and reduced R square.
				RsqChange <- RsqFull - RsqReduced
				res <- c(RsqReduced, RsqChange)
				names(res) <- c(paste("Rsquare Without ", paste0(x, collapse = " "), collapse = ""), "RsquareChange")
            }            
					
		}else{
            
            	#If effN == FALSE, we proceed with the overall n.
            
           #Adjusted R-square calculations.
            	multiplierFull <- (n-1)/(n - pFull - 1)
			multiplierRed <- (n-1)/(n - pRed - 1)
			RsqReduced <- 1 - multiplierRed*(1 - RsqReduced)
			RsqFull <- 1 - multiplierFull*(1 - RsqFull)
			
			#R square change is difference between overall R square and reduced R square.
			RsqChange <- RsqFull - RsqReduced

			res <- c(RsqReduced, RsqChange)

			names(res) <- c(paste("Rsquare Without ", paste0(x, collapse = " "), collapse = ""), "RsquareChange")
		}
 
	}else{
		
		#Otherwise, simply calculate R-square change without the adjustment terms.
		RsqChange <- RsqFull - RsqReduced
		res <- c(RsqReduced, RsqChange)
		names(res) <- c(paste("Rsquare Without ", paste0(x, collapse = " "), collapse = ""), "RsquareChange")
	}
	
	#if silent printing is not requested, print the result.
	if(!silent) print(round(res, 2))
	
	#And return the object.
	invisible(res)
}

rsquareCalcMG <- function(model, y, x){
	#model is a model fit by lavaan
	#y is a character string specifying the name of structural outcome of interest.
	#x is a vector of one or more character strings specifying the name(s) of the target predictor(s) of interest, to be omitted from the reduced model.
	
	require(lavaan)
	#parameter estimates
	pe <- parameterEstimates(model, standardized = TRUE, rsquare = TRUE)
	
	#rsqc = modified rsquareCalc() function, nested within this function.
	rsqc <- function(peObj, mod){
		
		#correlation matrix of all variables
		RmatTest <- lavInspect(mod, what = "cor.all")
		
		#Is RmatTest a list of correlation matrices for each of G groups?
		if(is.list(RmatTest)){
			
			#If so, this will be applied to each unique block using lapply() below.
			#Retrieve the correlation matrix corresponding to the unique block.
			Rmat <- lavInspect(mod, what = "cor.all")[[peObj$block[1]]]
			
		}else{
			
			#Otherwise, there is only one matrix.
			Rmat <- lavInspect(mod, what = "cor.all")
			
		}
		
		#regression coefficients
		Beta <- peObj[peObj$lhs == y & peObj$op == "~", ]
		
		#names of X variables NOT specified in x (reduced model x's).
		otherXnames <- Beta[!(Beta$rhs %in% x), "rhs"]
				
		#Grab correlation matrix of other Xs.
		Rxx <- Rmat[otherXnames, otherXnames, drop = FALSE]
		
		#Inverse X cor mat.
		RxxInv <- solve(Rxx)
		
		#vector of xy correlations.
		Rxy <- Rmat[otherXnames,y, drop = FALSE] #this way preserves the correct order of the other x names
		
		#compute new betas as they would have been without the variables in x included in the model.
		#beta = RxxInv%*Rxy
		BetaNew <- RxxInv%*%Rxy
		
		#R square of the submodel
		Rsq <- t(BetaNew)%*%Rxy
		
		#R square change is difference between overall R square and reduced model R Square.
		RsqChange <- peObj[peObj$lhs == y & peObj$op == "r2", "est"] - Rsq
		
		if(!is.null(peObj$block)){
			
			#If the analysis is multiple group, add group name as part of the result vector.
			res <- c(peObj$block[1], Rsq, RsqChange)
			names(res) <- c(paste0("Group", peObj$block[1]), "RsquareReduced", "RsquareChange")
			
		}else{
			
			#Otherwise, just include Reduced R-squared and R-square Change.
			res <- c(Rsq, RsqChange)
			names(res) <- c(paste("Rsquare Without ", paste0(x, collapse = " "), collapse = " "), "RsquareChange")
			
		}
		return(res)
	}
	
	#Is analysis MG?
	if(is.null(pe$block)){
		
		#If not, just calculate R-squared change for the single model.
		Result <- rsqc(pe, model)
		
		#Results for pretty printing (rounded to 2 decimals).
		resPrint <- Result
		resPrint <- round(resPrint, 2)
		
		#Rename Result columns.
		names(Result) <- c("RsquareReduced", "RsquareChange")
		
	}else{
		
		#If MG, then ...
		#Split the pe object into groups (blocks).
		peList <- split(pe, pe$block)
		#Calculate Change in R-squared for each group (block).
		ResultList <- lapply(peList, rsqc, model)
		#Reformulate the list as a table.
		
		for(i in 1:length(ResultList)){
			if(i == 1) Result <- ResultList[[i]] else Result <- rbind(Result, ResultList[[i]])
		}
		
		#Set group colnames as "Group".
		colnames(Result)[1] <- "Group"
		#Ensure the table is saved as a data frame rather than a matrix.
		Result <- data.frame(Result)
		#Add group labels, for convenience.
		Result$GroupLabel <- model@Data@group.label
		#Change rownames to numbers, per usual.
		rownames(Result) <- as.character(1:nrow(Result))
		#Reorder columns
		Result <- Result[,c(1, 4, 2:3)]
		
		#Results for pretty printing (numeric columns rounded to 2 decimals).
		resPrint <- Result
		resPrint[,3:4] <- round(resPrint[,3:4], 2)
	}
	
	
	#Print pretty version.
	print(resPrint)
	
	#Return full version (not rounded to 2 decimals).
	invisible(Result)
}


rsquareCalc.Boot <- function(modelSyn, Data, Y, X, nboot = 1000, FUN = lavaan::sem, miss = FALSE, adj = FALSE, effN = FALSE, postcheck = TRUE, seed = NULL, parallel = FALSE, ncpus = 2, conflevel = .95, bc = FALSE, ...){
	#modelSyn is lavaan model syntax.
	#Data is the dataset to be used in fitting the model and bootstrapping.
	#Y is a character vector of length 1 containing the name of the outcome (Y) variable in the structural regression of interest.
	#X is a character vector of length >= 1 indicating the target predictor(s) omitted from the reduced model.
	#nboot is the number of bootstrap resamples. Defaults to 1000, same as lavaan's bootstrapping default.
	#FUN is the function from the lavaan package to be used in fitting the model. sem() is the default.
	#miss: is there missing data and, if so is fiml estimation desired, with fixed.x set to TRUE?
	#adj: do you want to calculate adjusted rather than unadjusted R2 and R2 change? Defaults to FALSE.
	#effN: if TRUE, N in the adjusted R-square calculation is set to the lowest effective N in the structural regression. Defaults to FALSE.
	#postcheck: if TRUE (default), bootstrap estimates are discarded from analyses of bootstrapped datasets that produce inadmissible solutions/Heywood cases.
	#seed: defaults to NULL. If a seed is specified, the pseudo random number generating seed is set using set.seed(seed).
	#parallel: is parallel processing desired? Defaults to FALSE.
	#ncpus: if parallel = TRUE, how many cores do you wish to call for parallel processing? Defaults to 2.
	#conflevel: confidence level for the interval. Defaults to .95 for a 95% bootstrap CI.
	#bc: bias-corrected boostrap confidence interval is requested if bc = TRUE. Default is FALSE, yielding a percentile bootstrap CI.
	#... optional arguments to FUN.
	
	#require the lavaan package
	require("lavaan")
	if(!is.character(Y)|length(Y) != 1) stop("Y must be a character vector of length 1, specifying the name of the DV in the (manifest or latent variable) regression of interest!")
	
	#N is the number of rows in the dataset.
	N <- nrow(Data)
	
	#If a user-specified seed is provided, set the seed to ensure a replicable analysis.
	if(!is.null(seed)) set.seed(seed)
	
	#Matrix of bootstrap row indices.
	#nrow = N, and every column contains the indices for one bootstrap resample.
	#Pre-generating the indices in this way ensures that the random number seed
	#provides replicable results, even when using parallel processing.
	indMat <- matrix(sample(1:N, nboot*N, TRUE), nrow = N, ncol = nboot)
	
	#Function to pass to sapply() or sfClusterApplyLB():
	r2bootcalc <- function(colnum){
		#colnum is an index giving the column of the indMat being used
		#at a particular iteration.
		
		#Print colnum, giving the bootstrap sample currently under analysis.
		print(paste0("Bootstrap Sample: ", colnum))
		
		#Take column of indices:
		ind <- indMat[,colnum]
		
		#Subset data:
		bootData <- Data[ind,]
		
		if(miss == TRUE){
			
			#If missing data estimation is requested, set missing = "fiml" and fixed.x = FALSE:
			bootMod <- try(FUN(model = modelSyn, data = bootData, missing = "fiml", fixed.x = FALSE, ...))
			
		}else{
				
			#Try to fit user model to bootstrapped data:
			bootMod <- try(FUN(model = modelSyn, data = bootData, ...))
				
		}
		
		
		if(class(bootMod) == "try-error" | lavInspect(bootMod, what = "converged") == FALSE){
			
			#If the model does not converge, results is NA:
			R2ch <- NA
			
		}else if(postcheck){
			
			if(!lavInspect(bootMod, what = "post.check")){
				#If the user requests postcheck option (the default), models with Heywood cases are discarded:
				R2ch <- NA
				
			}else{
					
					#Generate rsquare change:
					R2ch <- rsquareCalc(bootMod, Y, X, adj, effN, silent = TRUE)["RsquareChange"]
			}
				
		}else{
			
			#In all other conditions, generate rsquare change:
			R2ch <- rsquareCalc(bootMod, Y, X, adj, effN, silent = TRUE)["RsquareChange"]
			
		}
		
		#Return R2ch:
		return(R2ch)
	}
	
	
	if(parallel){
		#If parallel computing is required, load the snowfall package:
		require("snowfall")
		
		#Initialize parallel computing:
		sfInit(parallel = TRUE, cpus = ncpus)
		
		#Load lavaan on the nodes:
		sfLibrary(lavaan)
		
		#Export relevant objects to the nodes:
		sfExport(list = c("modelSyn", "Data", "Y", "X", "nboot", "FUN", "miss", "adj", "effN", "postcheck", "indMat", "rsquareCalc"))
		
		#deltaR2Boot is the vector resulting from unlisting the cluster lapply (load balanced) function results:
		deltaR2Boot <- unlist(sfClusterApplyLB(x = 1:nboot, fun = r2bootcalc))
		
		#Stop parallel computing
		sfStop()
	}else{
		
		#For non-parallel computing, simply use sapply:
		deltaR2Boot <- sapply(X = 1:nboot, FUN = r2bootcalc)
		
	}
	
	#Internal function from the boot package:
	norm.inter <- function (t, alpha){
	    t <- t[is.finite(t)]
	    R <- length(t)
	    rk <- (R + 1) * alpha
	    if (!all(rk > 1 & rk < R)) 
	        warning("extreme order statistics used as endpoints")
	    k <- trunc(rk)
	    inds <- seq_along(k)
	    out <- inds
	    kvs <- k[k > 0 & k < R]
	    tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs + 
	        1))))
	    ints <- (k == rk)
	    if (any(ints)) 
	        out[inds[ints]] <- tstar[k[inds[ints]]]
	    out[k == 0] <- tstar[1L]
	    out[k == R] <- tstar[R]
	    not <- function(v) xor(rep(TRUE, length(v)), v)
	    temp <- inds[not(ints) & k != 0 & k != R]
	    temp1 <- qnorm(alpha[temp])
	    temp2 <- qnorm(k[temp]/(R + 1))
	    temp3 <- qnorm((k[temp] + 1)/(R + 1))
	    tk <- tstar[k[temp]]
	    tk1 <- tstar[k[temp] + 1L]
	    out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 - 
	        tk)
	    cbind(round(rk, 2), out)
	}
	
	#Fit model on original sample data with missing data specifications, as above:
	if(miss){
		estMod <- try(FUN(model = modelSyn, data = Data, missing = "fiml", fixed.x = FALSE, ...))
	}else{
		estMod <- try(FUN(model = modelSyn, data = Data, ...))
	}
	
	#r2change est from original sample data:
	Est <- rsquareCalc(estMod, Y, X, adj, effN, silent = TRUE)["RsquareChange"]
	
	#Bootstrap standard error is standard deviation of bootstrap estimates:
	bootSE <- sd(deltaR2Boot, na.rm = TRUE)
	
	###Next lines are adapted from parameterEstimates() function in lavaan:
	alpha <- (1 + c(-conflevel, conflevel))/2
	if(bc){
		zalpha <- qnorm(alpha)
		w <- qnorm(sum(na.omit(deltaR2Boot) < Est)/length(na.omit(deltaR2Boot)))
		a <- 0
		adj.alpha <- pnorm(w + (w + zalpha)/(1 - a *(w + zalpha)))
		qq <- norm.inter(na.omit(deltaR2Boot), adj.alpha)
		
	}else{
		qq <- norm.inter(na.omit(deltaR2Boot), alpha)
	}
	
	boot.CI.Lower = qq[1,2]
	boot.CI.Upper = qq[2,2]
	###
	
	names(boot.CI.Lower) <- names(boot.CI.Upper) <- NULL
	
	#Number of converged bootstrap replications = number of elements of bootstrap result vector that are not NA:
	converged <- length(which(!is.na(deltaR2Boot)))
	
	#save results:
	res <- data.frame(Est = Est, bootSE = bootSE, boot.CI.Lower = boot.CI.Lower, boot.CI.Upper = boot.CI.Upper, bootConvergence = converged)
	rownames(res) <- NULL
	#return results:
	return(res)
}
