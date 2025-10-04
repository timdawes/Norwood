# "Predicting clinical trajectory after the Norwood procedure: focus on aortic and pulmonary re-interventions"
# Authors: Haapanen H, Dawes TJW, Brown K, Giardini A, Dedieu N, Shetty P, Tsang V, Kostolny M.


# Functions



# Load packages

options(warn=-1)

suppressPackageStartupMessages({
  library(readxl)
  library(survival)
  library(Hmisc)
  library(lubridate)
  library(survminer)
  library(ggplot2)
  library(km.ci)
  library(ggsurvfit)
  library(mstate)
  library(flexsurv)
  library(bcv)
  library(dplyr)
  library(reshape)
  library(tidyverse)
  library(ggstream)
  library(missForest)
  library(ggpubr)
  library(mgcv)
  library(forecast)
  library(fabricatr)
  library(data.table)
  library(mltools)
})


rescale.yaxis<- function(x) {ifelse(x<1, x, 1+(x-1)/10)}

smoother<- function(new.rows, bw)
      {
        new.rows$mean_smoothed<- ksmooth(new.rows$time, new.rows$mean, "normal",bandwidth = bw)$y
        new.rows$ymin_smoothed<- ksmooth(new.rows$time, new.rows$ymin, "normal",bandwidth = bw)$y
        new.rows$ymax_smoothed<- ksmooth(new.rows$time, new.rows$ymax, "normal",bandwidth = bw)$y
        
        new.rows$ymin_smoothed[new.rows$ymin_smoothed < 0]<- 0
        new.rows$ymax_smoothed[new.rows$ymax_smoothed > 100]<- 100
        
        return(new.rows)  
      }


mean.sd.after.boxcox<- function(x, lambda)
{
  boxcox.transformed<- (x^lambda - 1) / lambda
  mean.bc<- mean(boxcox.transformed)
  upperCI.bc<- mean.bc + 1.96*sd(boxcox.transformed)
  lowerCI.bc<- mean.bc - 1.96*sd(boxcox.transformed)
  
  mean.inv.boxcox.transformed<- (mean.bc * lambda + 1) ^ (1/lambda)
  upperCI.inv.boxcox.transformed<- (upperCI.bc * lambda + 1) ^ (1/lambda)
  lowerCI.inv.boxcox.transformed<- (lowerCI.bc * lambda + 1) ^ (1/lambda)
  
  return (c(mean.inv.boxcox.transformed, lowerCI.inv.boxcox.transformed, upperCI.inv.boxcox.transformed))
}



expand.covs.msdata.TD <- function(data, covs, append=TRUE, longnames=TRUE, ...)
{
  if (!inherits(data, "msdata"))
    stop("'data' must be an 'msdata' object")
  trans <- attr(data, "trans")
  data <- as.data.frame(data)
  trans2 <- to.trans2(trans)
  K <- nrow(trans2)
  if (is.character(covs)) form1 <- as.formula(paste("~ ",paste(covs,collapse=" + "))) else form1 <- as.formula(covs)
  # going to apply model.matrix, but NA's are not allowed, so have to deal with that
  mm4 <- NULL
  for (j in 1:length(covs)) {
    wh <- which(!is.na(data[[covs[j]]]))
    form1 <- as.formula(paste("~ ",covs[j]))
    mm <- model.matrix(form1,data=data)
    mm <- data.frame(mm)
    mm <- mm[,-1,drop=FALSE]
    if (!longnames) {
      nc <- ncol(mm)
      if (nc==1) names(mm) <- covs[j]
      else names(mm) <- paste(covs[j],1:nc,sep="")
    }
    nms <- names(mm)
    ms <- data.frame(trans=data[["trans"]])
    ms$trans <- factor(ms$trans, levels = as.character(1:K)) # This is the line which differs from the original function in stating: 'levels = as.character(1:K)'
    ms <- cbind(ms[wh,,drop=FALSE],mm)
    mm2 <- model.matrix(as.formula(paste("~ (",paste(nms,collapse=" + "),"):trans")),data=ms)[,-1]
    mm3 <- matrix(NA,nrow(data),ncol(mm2))
    mm3[wh,] <- mm2
    mm3 <- data.frame(mm3)
    nms <- as.vector(t(outer(nms,1:K,"paste",sep=".")))
    names(mm3) <- nms
    if (j==1) mm4 <- mm3 else mm4 <- cbind(mm4,mm3)
  }
  
  if (!append) return(mm4) else {
    if (!all(is.na(match(names(data),nms))))
      warning("One or more names of appended data already in data!")
    mm4 <- cbind(data,mm4)
  }
  
  class(mm4) <- c("msdata", "data.frame")
  attr(mm4, "trans")<- trans
  return(mm4)
}





draw.survival.groups<- function (progression.mat)
        {
        # Create stacked data for areas
        stacked_values <- t(apply(progression.mat[,-1], 1, cumsum))
        times<- progression.mat[,1]
        
        
        # Base plot
        plot(NA, xlab = "", ylab = "", xlim=range(times), ylim=c(0, max(stacked_values)), type="n", axes=FALSE)
        axis(side = 1, at = seq(0,20,5), labels = seq(0,20,5), line = 0.5, mgp = c(3,0.5,0), las=1)
        axis(side = 2, at = seq(0,1,0.2), labels = seq(0,1,0.2), line = 0.5, mgp = c(3,0.5,0), las=2)
        
        mtext("Time after Norwood (years)", side = 1, line = 2)
        mtext("Probability", side = 2, line = 2.5)
        
        
        # Add stacked areas
        for (i in 1:4) {
          if (i == 1) {polygon(c(times, rev(times)), c(rep(0, length(times)), rev(stacked_values[, i])), col=rainbow(20)[c(1,4,9,14)][i], border=NA)} else 
          {polygon(c(times, rev(times)), c(stacked_values[, i - 1], rev(stacked_values[, i])), col=rainbow(20)[c(1,4,9,14)][i], border=NA)}
        }
        
        # Add lines on top of the stacked areas
        matlines(times, stacked_values, lty=1, lwd=1.5, col="black")
        
        stacked.lines<- c(0, stacked_values[which(times == 10),])
        label.heights<- 0.5*(stacked.lines[-1] + stacked.lines[-5])
        sapply(1:4, function(x) {text(10, label.heights[x], c("Direct pathway", "Directplus pathway", "Indirect pathway", "Dead or transplanted")[x])})
}




redrank.TD <- function(redrank, full = ~1, data, R,
                       strata = NULL, Gamma.start, method = "breslow", eps = 1e-5, 
                       print.level = 1, max.iter = 100, remove.NAs = FALSE, scale = TRUE)
        {
                  if (!inherits(data, "msdata"))
                    stop("'data' must be an 'msdata' object")
                  
                  par(mfrow=c(1,1))
                  trans <- attr(data, "trans")
                  
                  # Now we need only the data
                        data <- as.data.frame(data)
                  
                  # Set variables as categorical or continuous/ordinal accordingly
                        cols.factors<- c("DGcat", "SysV", "AVVRcat", "SanovsBT", "Interdigitation", "N1ECMO")
                        for (i in cols.factors) {data[,i]<- factor(data[,i], levels = c(0,1))}
                        
                        cols.continuous<- c("Weight","SystV_f")
                        for (i in cols.continuous) {data[,i]<- as.numeric(data[,i])}
                        
                        
                  # Get model matrices of reduced rank and full rank parts
                        mmrr <- model.matrix(redrank, data=data)
                        mmrr <- mmrr[,-1,drop=FALSE] # without intercept
                        
                  # Scale the covariate data to improve fitting
                        if (scale == TRUE) {
                            mmrr<- scale(mmrr)
                            mmrr.center<- attr(mmrr, "scaled:center")
                            mmrr.scale<- attr(mmrr, "scaled:scale")} else {mmrr.center<- mmrr.scale<- NA}
                            
                        p <- ncol(mmrr)
                        if (p==0) stop("Empty reduced rank part; please consider full model")
                        mmrr <- data.frame(mmrr)
                        covs <- names(mmrr)
                        mmfull <- model.matrix(full, data=data)
                        mmfull <- mmfull[,-1,drop=FALSE] # without intercept
                        p2 <- ncol(mmfull)
                        
                  # Construct working data
                        if (p2>0) {
                          mmfull <- data.frame(mmfull)
                          fullcovs <- names(mmfull)
                          rrdata <- as.data.frame(data[,c("id","from","to","trans","Tstart","Tstop","time","status")])
                          rrdata <- cbind(rrdata,mmrr,mmfull)
                          cols <- 8 + (1:p) # Finds which columns correspond to the covariates (Z)
                          fullcols <- 8 + p + (1:p2)
                        } else {
                          rrdata <- as.data.frame(data[,c("id","from","to","trans","Tstart","Tstop","time","status")])
                          rrdata <- cbind(rrdata,mmrr)
                          cols <- 8 + (1:p)
                          fullcols <- NULL        
                        }
                        
                  # Get and store whether clock is forward or reset
                        cx <- coxph(redrank, data=data)
                        if (attr(cx$y, "type") == "counting") clock <- "forward" else if (attr(cx$y, "type") == "right") clock <- "reset" else stop("Surv object should be either of type 'counting' or 'right'")
                        
                  # Preparations for iterative algorithm
                        trans2 <- to.trans2(trans)
                        K <- nrow(trans2) # K = number of transitions
                        if (!is.null(dimnames(trans)))
                          tnames <- paste(trans2$fromname,"->",trans2$toname) else tnames <- as.character(1:K)
                        
                        ### add to the data set R replicates of columns with covariates Z_1...Z_p
                        
                        colsR <- matrix(0,R,p) # colsR is a matrix containing reduced-rank covariates
                        
                        for (r in 1:R) {
                          ncd <- ncol(rrdata)
                          rrdata <- cbind(rrdata,rrdata[,cols])
                          colsR[r,] <- ((ncd+1):(ncd+p))
                          names(rrdata)[((ncd+1):(ncd+p))] <- paste(covs, as.character(r), sep=".rr")
                        }
                        
                        Delta<- 100
                        old.gammas<- matrix(0, nrow=K, ncol=1)
                        old.alphas<- matrix(0, nrow=p, ncol=1)
                        par(mfrow=c(1,3))
                        
                        while(abs(Delta) > eps) {
                          
                            Gamma.start <- matrix(rnorm(R*K, mean=0, sd=0.1),R,K) # Gamma = influence of the coefficients on each transition (dimensions: R x K)
                            Gamma.iter <- Gamma.start
                            iter <- 1
                            prev.loglik <- 0
                            loglik <- 100
                            old.deltas<- 1
                            Delta <- loglik - prev.loglik
                            
                        
                                  while(abs(Delta) > eps && iter < max.iter) {
                                    if (print.level > 0) {
                                      cat("\n\nIteration", iter, "... \n")
                                      flush.console()
                                    }
                                    iter <- iter + 1
                                    class(rrdata) <- c("msdata", "data.frame")
                                    attr(rrdata, "trans")<- trans
                                
                                    ms.it <- redrank.iter(data = rrdata, cols = cols, colsR = colsR, fullcols = fullcols, R = R, clock = clock, strata = strata, Gamma.iter = Gamma.iter, print.level = print.level - 1)
                                    
                                    Gamma.iter <- ms.it$Gamma
                                    if (print.level > 0) cat("\nGamma = ", round(Gamma.iter, 5))
                                    prev.loglik<- loglik
                                    loglik<- ms.it$loglik
                                    Delta <-  - loglik + prev.loglik
                                    
                                    old.deltas<- c(old.deltas, Delta)
                                    old.gammas<- cbind(old.gammas, matrix(Gamma.iter, ncol=1))
                                    old.alphas<- cbind(old.alphas, matrix(ms.it$Alpha, ncol=1))
                                    
                                    
                                    if (print.level > 0) {
                                      cat("\nPrevious loglik =", round(prev.loglik, 5), ", present loglik =", round(loglik, 5),
                                          " Delta = ", round(Delta, 8))
                                      plot(1:iter, -log10(abs(old.deltas)), type='l', lwd=4, col="blue", ylim=c(min(-log10(abs(old.deltas))),-log10(eps)+0.5))
                                      matplot(1:iter, t(old.gammas), type='l', lwd=3)
                                      matplot(1:iter, t(old.alphas), type='l', lwd=3)
                                      abline(h=5, lty=2, col="red", lwd=3)
                                      flush.console()
                                    }
                                    
                                    # Rescale Gamma.iter
                                                Alpha <- ms.it$Alpha
                                                B <- Alpha %*% ms.it$Gamma
                                                svd.B <- svd(impute.svd(B)$x) # Changed from svd() to impute.svd() by TD
                                                l<- length(svd.B$d)
                                                
                                            # Reconstruct Gamma
                                                Gamma.iter <- (diag(sqrt(svd.B$d), nrow=l, ncol=l) %*% t(svd.B$v))[1:R,]
                                                Gamma.iter <- matrix(Gamma.iter,R,K)
                                                rnames <- paste("r",1:R,sep="")
                                                dimnames(Gamma.iter) <- list(rnames,tnames)
                                                
                                            # Reconstruct Alpha
                                                Alpha <- (svd.B$u %*% diag(sqrt(svd.B$d), nrow=l, ncol=l))[,1:R]
                                            
                                            # Normalise Alpha to have unit-length columns
                                                if (R>1) norm.Alpha <- apply(Alpha^2,2, function(x) sqrt(sum(x))) else norm.Alpha <- sqrt(sum(Alpha^2))
                                                norm.Alpha.mat <- matrix(norm.Alpha,p,R,byrow=TRUE)
                                                Alpha <- Alpha/norm.Alpha.mat
                                                dimnames(Alpha) <- list(covs,rnames)
                                                AlphaX <-  as.matrix(rrdata[,cols]) %*% Alpha
                                                
                                            # Gamma re-scaled to maintain the scale of B
                                                Gamma.iter <- Gamma.iter * matrix(norm.Alpha,R,K)
                                    
                                    
                                    
                                    
                                  }
                        }
                
                
                
         
                svd.B <- svd(impute.svd(B)$x) 
                l<- length(svd.B$d)
                
              # Reconstruct Gamma
                      Gamma.final <- (diag(sqrt(svd.B$d), nrow=l, ncol=l) %*% t(svd.B$v))[1:R,]
                      Gamma.final <- matrix(Gamma.final,R,K)
                      rnames <- paste("r",1:R,sep="")
                      dimnames(Gamma.final) <- list(rnames,tnames)
                      
              # Reconstruct Alpha
                      Alpha <- (svd.B$u %*% diag(sqrt(svd.B$d), nrow=l, ncol=l))[,1:R]
                      
                      # Normalise Alpha to have unit-length columns
                            if (R>1) norm.Alpha <- apply(Alpha^2,2, function(x) sqrt(sum(x))) else norm.Alpha <- sqrt(sum(Alpha^2))
                            norm.Alpha.mat <- matrix(norm.Alpha,p,R,byrow=TRUE)
                            Alpha <- Alpha/norm.Alpha.mat
                            dimnames(Alpha) <- list(covs,rnames)
                            AlphaX <-  as.matrix(rrdata[,cols]) %*% Alpha
                            
              # Gamma re-scaled to maintain the scale of B
                      Gamma.final <- Gamma.final * matrix(norm.Alpha,R,K)
                      
                      
                dimnames(B) <- list(covs,tnames)
          return(list(Alpha = Alpha, Gamma = Gamma.final, Beta = B, Beta2 = ms.it$Beta2, cox.itr1 = ms.it$cox.itr1, cox.itr2 = ms.it$cox.itr2,
                      AlphaX = AlphaX, niter = iter, df = R*(p+K-R), loglik = loglik, negative.log10.pvals = ms.it$negative.log10.pvals,
                      scale = mmrr.scale, center = mmrr.center))
        }
        




redrank.iter <- function(data, cols, colsR, fullcols, R,
                         clock, strata, Gamma.iter, print.level = 0)
          {
            if (!inherits(data, "msdata"))
              stop("'data' must be an 'msdata' object")
            trans <- attr(data, "trans")
            data[,c(cols,colsR,fullcols)] <- unlist(data[,c(cols,colsR,fullcols)])
            nfull <- length(fullcols)
            covariates <- names(data)[cols]
            fullcovariates <- names(data)[fullcols]
            n <- nrow(data)
            K <- max(trans,na.rm=TRUE)
            p <- length(cols)
            
            # Scales the covariates using Gamma
            # This section adjusts the covariates for each transition 'k', applying the current estimate of 'Gamma' to scale covariates relevant to rank R
                  for (k in 1:K) { ### W[k,r] = Gamma.iter[r,k] * Z
                    wh <- which(data$trans == k)
                    for (r in 1:R) {
                      data[wh, colsR[r,]] <- Gamma.iter[r,k] * data[wh, colsR[r,]]
                    }
                  }
            data$trans <- factor(data$trans, levels = as.character(1:K)) ### lost in the unlist above
            covs.R <- names(data)[as.vector(t(colsR))]
            if (is.null(strata)) strata <- "trans"
            strata.expr <- paste("strata(",strata,")+")
            
            # First iteration estimates alpha by fitting a Cox model to scaled covariates stratified by transition
                  # Using ridge regression
                        if (clock=="forward")
                          expr1 <- paste("cox.itr1 <- coxph(Surv(Tstart,Tstop,status)~",
                                         strata.expr,
                                         "ridge(",
                                         paste(covs.R, collapse = ","), ", theta = 1.0)", sep="") else
                          expr1 <- paste("cox.itr1 <- coxph(Surv(time,status)~",
                                         strata.expr,
                                         paste(covs.R, collapse = "+"), sep="")
                
                        if (!is.null(fullcols)) expr1 <- paste(expr1, "+", paste(fullcovariates, collapse = "+"),sep="")
                        expr1 <- paste(expr1,", data=data, na.action=\"na.exclude\", )", sep = "")
                
               
                      cox.itr1 <- eval(parse(text = expr1, n = 1))
                      if (print.level > 0) {
                        cat("\nStratified Cox regression:\n\n")
                        print(cox.itr1)
                      }
                      loglik <- cox.itr1$loglik[2]
                      if (print.level > -1) cat("\nAlpha =", round(cox.itr1$coef, 5))
                    
                  # Alpha is the coefficients for the covariate vector Z - i.e. influence of each covariate on the reduced rank coefficient  
                      Alpha <- matrix(cox.itr1$coef[1:(p*R)],p,R) 
            
              # Second iteration estimates Gamma by constructing prognostic scores (Alpha . covariates)
                      
                      Beta2 <- cox.itr1$coef[-(1:(p*R))]
                      ncd <- ncol(data)
                      
                      # Multiply the covariates (Z) by Alpha for each patient and make a new column for each rank
                            data <- cbind(data,as.matrix(data[, cols]) %*% Alpha)
                            AlphaX.R <- paste("AlphaX",as.character(1:R),sep="") # Name it AlphaX.R (1->R)
                            names(data)[((ncd+1):(ncd+R))] <- AlphaX.R
                            
                      # Now expand the covariates for this covariate
                            attr(data, "trans") <- trans
                            class(data) <- c("msdata", "data.frame")
                            data <- expand.covs.msdata.TD(data,AlphaX.R)
                            AlphaX.RK <- names(data)[((ncd+R+1):(ncd+R+R*K))] 
                      
                      # Do ridge regression
                      if (clock=="forward")
                        expr2 <- paste("cox.itr2 <- coxph(Surv(Tstart,Tstop,status)~",
                                       strata.expr,
                                       "ridge(",
                                       paste(AlphaX.RK, collapse = ","), ", theta = 1.0)", sep="") else
                                         expr2 <- paste("cox.itr2 <- coxph(Surv(time,status)~", strata.expr, paste(AlphaX.RK, collapse = "+"),sep="")
                    
                       
                      if (!is.null(fullcols)) expr2 <- paste(expr2, "+", paste(fullcovariates, collapse = "+"),sep="")
                      expr2 <- paste(expr2,", data=data, na.action=\"na.exclude\")", sep = "")
                      cox.itr2 <- eval(parse(text = expr2, n = 1))
                      
            # TD added
                  negative.log10.pvals<- -log10(coefficients(summary(cox.itr2))[,'p'])
                  negative.log10.pvals[is.na(negative.log10.pvals)]<- 0
                  
            if (print.level > 0) {
              cat("\n\nCox regression on scores\n\n")
              print(cox.itr2)
            }  
            
            # Extract the Gamma coefficients for each transition
                  Gamma.iter <- t(matrix(cox.itr2$coef[1:(K*R)],K,R))
                  Beta2 <- cox.itr2$coef[-(1:(K*R))]
                  
            return(list(Gamma = Gamma.iter, Alpha = Alpha, Beta2 = Beta2, loglik = loglik, cox.itr1 = cox.itr1, cox.itr2 = cox.itr2, negative.log10.pvals = negative.log10.pvals))
          }




mssample.TD <- function(Haz, trans, history=list(state=1,time=0,tstate=NULL),
                       beta.state=NULL,
                       clock=c("forward","reset"),
                       output=c("state","path","data"), tvec, cens=NULL, M=10, do.trace=NULL)
{ 

  
  output <- match.arg(output)
  clock <- match.arg(clock)
  K <- dim(trans)[1]
  trans2 <- to.trans2(trans)
  ntrans <- nrow(trans2)
  if (length(history$state)==1) history$state <- rep(history$state,M)
  if (length(history$time)==1) history$time <- rep(history$time,M)
  if (length(history$state)!=length(history$time)) stop("lengths of history$state and history$time differ")
  if (!is.null(history$tstate)) { # then should be either length K or dim 
    if (is.vector(history$tstate)) 
      if (length(history$tstate) != K) stop("length of history$tstate should equal no of states")
    else history$tstate <- matrix(history$tstate,K,M)
    if (is.null(beta.state)) stop("beta.state should be specified when history$tstate not null")
  }
  if (!is.null(beta.state))
    if (any(dim(beta.state) != c(K,ntrans)))
      stop("incorrect dimension of beta.state")
  if (output=="state") # to contain sum
    res <- matrix(0, length(tvec), K) else if (output=="path") {# to contain sum
    thepaths <- paths(trans)
    L <- nrow(thepaths)
    res <- matrix(0, length(tvec), L)
  } else res <- NULL
  for (m in 1:M) {
    if (!is.null(history$tstate))
      res1 <- mssample1.TD(Haz, trans, history=list(state=history$state[m], time=history$time[m], tstate=history$tstate[,m]), beta.state=beta.state, clock=clock, output=output, tvec=tvec, cens=cens) else
      res1 <- mssample1.TD(Haz, trans, history=list(state=history$state[m], time=history$time[m], tstate=rep(0,K)), beta.state=beta.state, clock=clock, output=output, tvec=tvec, cens=cens)
    if (output=="data") {
      res1[,1] <- m
      res <- rbind(res,res1)
    }
    else res <- res + res1
    if (!is.null(do.trace)) if (m %% do.trace == 0) {
      cat("Replication",m,"finished at",date(),"\n")
      flush.console()
    }
  }
  if (output=="state") {
    res <- data.frame(cbind(tvec,res/M))
    names(res) <- c("time",paste("pstate",1:K,sep=""))
  } else
  
  if (output=="path") {
    res <- data.frame(cbind(tvec,res/M))
    names(res) <- c("time",paste("ppath",1:L,sep=""))
  } else
  
  if (output=="data") {
    res <- data.frame(res)
    names(res) <- c("id","Tstart","Tstop","duration","from","to","status","trans")
    attr(res, "trans") <- trans
    class(res) <- c("msdata", "data.frame")
  }
  return(res)  
}



              
mssample1.TD <- function(Haz, trans, history, beta.state, clock, output, tvec, cens)
{ 
          if (!is.null(cens)) {
            pcens <- diff(c(0,1-cens$surv))
            idx <- sample(1:length(cens$time), size=1, prob=pcens)
            fut <- cens$time[idx]
            censtime <- list(time=fut, jump=ifelse(idx>1,cens$Haz[idx]-cens$Haz[idx-1],cens$Haz[idx]))
          } else censtime <- NULL
          K <- dim(trans)[1]
          trans2 <- to.trans2(trans)
          from <- to <- history$state
          tcond <- t0 <- Tstart <- history$time
          if (output=="state")
            res <- matrix(0, length(tvec), K) else if (output=="path") {
            thepaths <- paths(trans)
            path <- c(to, rep(NA,ncol(thepaths)-1))
            res <- matrix(0, length(tvec), nrow(thepaths))
          } else res <- NULL
  ### keep track of when states were visited
  tstates <- history$tstate
  while (!is.na(to)) {
    from <- to
    nstates <- trans[from,]
    transs <- nstates[!is.na(nstates)]
    allto <- which(!is.na(nstates))
    ntr <- length(transs)
    
    if (ntr!=0 ) { # if not yet in absorbing state
      transnos <- transs
      for (tr in 1:ntr)
        Haz$Haz[Haz$trans==transnos[tr]] <-
          exp(sum(beta.state[,transnos[tr]]*tstates)) *
          Haz$Haz[Haz$trans==transnos[tr]]
      whh <- which(!is.na(match(Haz$trans,transnos)))
      if (clock=="forward") {
        crs <- crsample.TD(Haz[whh,], tcond, censtime)
        tcond <- Tstop <- crs$t
      }
      else {
        crs <- crsample.TD(Haz[whh,], t0, censtime)
        t0 <- 0
        tcond <- Tstop <- crs$t + tcond
      }
      transno <- crs$trans
      if (is.na(transno)) to <- NA
      else {
        to <- trans2$to[transno]
        tstates[to] <- Tstop
      }
      if (output=="state") {
        res[((tvec>=Tstart)&(tvec<Tstop)),from] <- 1
        Tstart <- Tstop
      }
      else if (output=="path") {
        idx <- which(apply(thepaths,1,function(x) identical(x,path)))
        res[((tvec>=Tstart)&(tvec<Tstop)),idx] <- 1
        path[which(is.na(path))[1]] <- to
        Tstart <- Tstop
      }
      else {
        res1 <- matrix(c(rep(NA,ntr),rep(Tstart,ntr),rep(Tstop,ntr),rep(Tstop-Tstart,ntr),rep(from,ntr),allto,rep(0,2*ntr)),ntr,8)
        res1[res1[,6]==to,7] <- 1
        res1[,8] <- trans[from,allto] # trans
        Tstart <- Tstop
        res <- rbind(res,res1)
      }
    } 
    else { # in absorbing state:
      to <- NA
      if (output=="state") {
        res[tvec>=Tstart,from] <- 1
      }
      else if (output=="path") {
        idx <- which(apply(thepaths,1,function(x) identical(x,path)))
        res[tvec>=Tstart,idx] <- 1
        path[which(is.na(path))[1]] <- to
      }
      else {
        res1 <- matrix(c(rep(NA,ntr),rep(Tstart,ntr),rep(Tstop,ntr),rep(Tstop-Tstart,ntr),rep(from,ntr),allto,rep(0,2*ntr)),ntr,8)
        res1[res1[,6]==to,7] <- 1
        res1[,8] <- trans[from,allto] # trans
        res <- rbind(res,res1)
      }
    }
  }
  return(res)  
}




crsample.TD <- function(Haz, tcond=0, censtime=NULL) ### censtime here is list with time and last jump (=pp.fut)
{
  if (is.null(censtime)) fut <- Inf else fut <- censtime$time
  transs <- Haz$trans
  transun <- unique(transs)
  K <- length(transun)
  tt <- sort(unique(Haz$time))
  n <- length(tt)
  cim <- matrix(NA, n, 3*K+4)
  ci <- as.data.frame(cim)
  names(ci)[1] <- "time"
  names(ci)[2:(K+1)] <- paste("Haz",as.character(1:K), sep="")
  names(ci)[(K+2):(2*K+1)] <- paste("haz",as.character(1:K), sep="")
  names(ci)[(2*K+2):(3*K+1)] <- paste("CI",as.character(1:K), sep="")
  names(ci)[3*K+2] <- "hazsum"
  names(ci)[3*K+3] <- "Hazsum"
  names(ci)[3*K+4] <- "S0"
  ci$time <- tt
  for (k in 1:K) { # select the elements for each transition in the block
    wh <- which(Haz$trans==transun[k])
    idx <- match(Haz$time[wh],tt)
    ci[,k+1][idx] <- Haz$Haz[wh]
    ci[,k+1] <- NAfix(ci[,k+1],subst=0)
    ci[,K+1+k] <- diff(c(0,ci[,k+1])) # the hazard for transition k
  }
  ## select only t > tcond, and adjust the cumulative hazards
  ci <- ci[ci$time>tcond,]
  n <- nrow(ci)
  for (k in 1:K) # each cumulative hazard is adjusted (cumulative sum of hazard)
    ci[,k+1] <- cumsum(ci[,K+1+k])
  ### compute baseline survival S0, need to sum all columns Hazards
  if (K==1) ci$hazsum <- ci[,3] else ci$hazsum <- apply(ci[,((K+2):(2*K+1))],1,sum)
  ci$S0 <- cumprod(1-ci$hazsum)
  ci$Hazsum <- -log(ci$S0)
  non.na.rows<- which(is.na(ci$Hazsum) == FALSE)
  ci<- ci[non.na.rows,]
  nci <- nrow(ci)
  ### Follow Dabrowska in first sampling from the sum of the cumulative hazards
  k <- NA
  tsample <- Hazsample.TD(data.frame(time=ci$time, Haz=ci$Hazsum))
  ### If follow-up time is less than tsample, no transition is made
  if (fut < tsample) crt <- fut else { # else decide between transitions
    crt <- tsample
    if (fut>tsample) # the size of the hazard jumps of the transitions
      # determine the probabilities of choosing the transitions
    {
      k <- sample(1:K,size=1,prob=ci[which(ci$time==tsample),(K+2):(2*K+1)])
    }
    else # i.e. if there is fut==tsample, then the jump in the censoring distribution
      # also plays in the lottery, unless both are Inf, then it doesn't matter
      if (crt!=Inf) {
        k <- sample(c(1:K,NA),size=1,
                    prob=c(ci[which(ci$time==tsample),(K+2):(2*K+1)],censtime$jump))
      }
  }
  if (!is.na(k)) trans <- unique(Haz$trans)[k] else trans <- NA
  return(list(t=crt,trans=trans))
}



Hazsample.TD <- function(Haz, size=1, replace=TRUE)
{
 
  p <- diff(c(0,1-exp(-Haz$Haz)))
  p <- c(p, exp(-Haz$Haz[nrow(Haz)])) # add probability of sampling time=Inf
  p[p<0]<- 0
  if (length(c(Haz$time, Inf)) == 1) {p<- 1}
  return(sample(c(Haz$time, Inf), size=size, prob=p, replace=replace))
}




`NAfix` <- function(x, subst=-Inf) {
 
  spec <- max(x[!is.na(x)])+1
  x <- c(spec,x)
  while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)]-1]
  x[x==spec] <- subst
  x <- x[-1]
  x
}


msboot.TD<- function (theta, data, cfull.Allc.boot, id = "id", B = 5, verbose = 0) 
      {
  
        par(mfrow=c(1,1))
        trans <- attr(data, "trans")
        ids <- unique(data[[id]])
        n <- length(ids)
        th <- theta(data, cfull.Allc.boot)
        res <- matrix(NA, length(th), B)
        colnames(res)<- paste("rep",1:B,sep="")
        rownames(res)<- names(th)
        
        for (b in 1:B) {
          if (verbose > 0) {
            r<- max(c(na.omit(res[!is.infinite(c(res))]),1))
            matplot(1:ncol(res), t(res[1:nrow(res),]), ylim=c(-r,r))
            flush.console()
          }
          bootdata <- NULL
         
          repeat {
                      bids <- sample(ids, replace = TRUE)
                      bidxs <- unlist(sapply(bids, function(x) which(x == data[[id]])))
                      bootdata <- data[bidxs, ]
                      print(events(bootdata))
                      cat("applying theta ...")
                      try(thstar <- theta(bootdata, cfull.Allc.boot))
                      if (class(thstar)!="try-error") {break}
                }
          
          res[, b] <- thstar
        }
        
        if (verbose) 
          cat("\n")
        return(res)
}





# This function generates bootstrap samples for the later analysis
          generate.bootstrap.sample.with.covariates<- function(era, all.data, replace) { 
            
                # Build the complete model
                      cfull.Allc.boot <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans) + Weight + Ndays + DGcat + SysV + SystV_f + AVVRcat + SanovsBT + Interdigitation + N1ECMO, data = all.data, method = "breslow")
                
                # Choose the rows for BOOTSTRASPPING
                      if (era == 0) {era.rows<- 1:210} else {era.rows<- which(data$era == era)}
                      era.data = all.data[era.rows,]
                      
                      trans<- attr(all.data, "trans")
                      newdata<- era.data[, c("Weight","Ndays","DGcat","SysV","SystV_f","AVVRcat","SanovsBT","Interdigitation","N1ECMO")]
                      newdata.typical<- matrix(0, nrow = 1, ncol = ncol(newdata), dimnames = list(1, colnames(newdata)))
                      
                      for (i in 1:ncol(newdata)) {
                        t<- sort(table(na.omit(newdata[,i])), decreasing=T)
                        if (length(t) > 5) {newdata.typical[1,i]<- sample(na.omit(newdata[,i]), 1)} else {newdata.typical[1,i]<- sample(as.numeric(names(t)),1)}
                      } 
                      
                      number.of.trans<- 39
                      newdata.multi<- data.frame(do.call(rbind, replicate(number.of.trans, newdata.typical, simplify=FALSE)))
                      newdata.multi$trans<- newdata.multi$strata<- 1:number.of.trans
                      attr(newdata.multi, "trans") <- trans
                      
                      msf.Allc.boot<- msfit(cfull.Allc.boot, newdata = newdata.multi, trans = tmat.Allc)
                      
                      return(msf.Allc.boot)
                    }
          
          generate.bootstrap.sample<- function(era, msm2.Allc.covs, replace) {
            
                    if (era == 1) {ids<- msm2.Allc.covs$id[which(msm2.Allc.covs$era == 1)]}
                    if (era == 2) {ids<- msm2.Allc.covs$id[which(msm2.Allc.covs$era == 2)]}
                    if (era == 3) {ids<- msm2.Allc.covs$id[which(msm2.Allc.covs$era == 3)]}
                    if (era == 0) {ids<- 1:210}
                    ids<- unique(ids)
                    
                    s<- sample(ids, length(ids), replace=replace)
                    boot.rows<- which(msm2.Allc.covs$id %in% s)
                    cfull.Allc.boot <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msm2.Allc.covs[boot.rows,], method = "breslow")
                    msf.Allc.boot<- msfit(cfull.Allc.boot, vartype = "greenwood", trans = tmat.Allc)
                    
                    return(msf.Allc.boot)
          }
                    




