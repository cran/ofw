#
# Copyright (C) 2007 Kim-Anh Lê Cao, Patrick Chabrier, INRA,
# French National Institute for Agricultural Research.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

"evaluateCART" <- function(obj, ...) UseMethod("evaluateCART")

"evaluateCART.learnCART"<-
  function(obj,
           maxvar=15,
           nvar=obj$nclass+1,
           ntreeTest=100,
	   weight=FALSE,
	   ...) {
    
    x=obj$x
    y=obj$y
    nsample=obj$nsample
    nvariable=obj$nvariable
    Bsample=obj$Bsample

   

    error.boot=vector(length=maxvar)
    vectWeight=vector(length=nrow(x))
    mat.pred.inbag=matrix(nrow=nrow(x), ncol=Bsample)
    mat.pred.test=matrix(nrow=nrow(x), ncol=Bsample)

    err.inbag=vector(length=Bsample)
    err.test=vector(length=Bsample)

    #for the weights if needed
    numWeight=matrix(nrow=nlevels(y), ncol=Bsample)
    classWeight=matrix(nrow=nlevels(y),ncol=Bsample)
    sampleWeight=matrix(nrow=nrow(x),ncol=Bsample)

    if(obj$weight==T  && weight==T){
    sampleWeight=obj$weightingOption$sampleWeight  #if obj has been weighted
    classWeight=obj$weightingOption$classWeight 
    }	

    #check input
    if(maxvar >=nvariable) stop("maxvar shoud not be greater than the number of variables")
     if(nvar >=nvariable) stop("nvar shoud not be greater than the number of variables")

    ##compute the weights (if needed) for each sample
    if(weight==T && obj$weight==F){
	for(boot in 1:Bsample){
		numWeight[,boot]=summary(y[obj$matTrain[,boot]])
		classWeight[,boot]=nsample/numWeight[,boot]
		classWeight[,boot]=classWeight[,boot]/sum(classWeight[,boot])
		for(n in 1:nlevels(y)){
		sampleWeight[which(as.integer(y[obj$matTrain[,boot]])==n),boot]= classWeight[n,boot]
	}
	}

    }
	
    
    
    ##BIG LOOP ON MAXvar
    cat("\n  Calculating e632+ for each of the ",maxvar," variables \n")
    for(ngene in 1:maxvar){
      cat("\n  variable :",ngene,"\n") 

      ##LOOP on the bootstrap samples
      for(boot in 1:Bsample){
        #cat("\n  boot :",boot,"\n") 

        
	proba=obj$matProb[,boot]
	train=obj$matTrain[,boot]
	test=setdiff(c(1:nsample),train) 
	xprim=x[train,]
	yprim=y[train]
    
    	nclass <- length(levels(y))
	
        
	##define test data
    	xtest=x[test,]
	ytest=y[test]
	ntest <- nrow(xtest)
        
	#if (weight==T){ vectWeight=sampleWeight[,boot]}
	if (weight==T){ vectWeight=classWeight[,boot]}	
        
                                        #choix fixe, voir plus tard majority vote wins here
	cutoff <- rep(1 / nclass, nclass)
                                        #no class weight here (a voir plus tard)
	classwt <- rep(1, nclass)
        
	nrnodes <- 2 * trunc(nrow(xprim) / min(summary(yprim))) + 1
        

	## Compiled code expects variables in rows and observations in columns.
	xprim <- t(xprim)
	storage.mode(xprim) <- "double"

	xtest <- t(xtest)
        storage.mode(xtest) <- "double"


	error.test <- double(ntreeTest)
	error.inbag <- double(ntreeTest)

	pred.test=as.integer(numeric(ntest))
	pred.inbag=as.integer(numeric(nsample))

	learnout <- .C("classLearn",
                    x = xprim,
                    xdim = as.integer(c(nvariable, nsample)),
                    y = as.integer(yprim),
                    nclass = as.integer(nclass),
                    ntree = as.integer(ntreeTest),
                    nvar = as.integer(nvar),
		    ngene = as.integer(ngene),
                    classwt = as.double(classwt),
                    cutoff = as.double(cutoff),
                    counttr = integer(nclass * nsample),
                    nrnodes = as.integer(nrnodes),
                    ndbigtree = integer(ntreeTest),
                    nodestatus = integer(nrnodes),
                    bestvar = integer(nrnodes),
                    treemap = integer(2 * nrnodes),
                    nodepred = integer(nrnodes),
                    xbestsplit = double(nrnodes),
                    xts = as.double(xtest),
                    clts = as.integer(ytest),
                    nts = as.integer(ntest),
                    countts = double(nclass * ntest),
                    outclts = pred.test,
		    outcl=pred.inbag,
                    errts = error.test,
		    errin = error.inbag,
		    proba=as.double(proba),
		    weight=as.integer(weight),
		    vectWeight=as.double(vectWeight),
                    DUP=FALSE,
                    PACKAGE="ofw")[-1]   


	##we suppose for inbag data that prediction is the same for identical cases
	
	mat.pred.inbag[train,boot]=learnout$outcl
	mat.pred.test[test,boot]=learnout$outclts

	err.test[boot]=mean(learnout$errts)
	err.inbag[boot]=mean(learnout$errin)

	#cat("\n  outcl :",errorout$outcl,"\n") 


	} #fin boot

#this e632+ code comes from the ipred package from Thorsten and from the varselRF package from Diaz

	one=mean(err.test, na.rm=T)
	#one <- mean(apply(cbind(mat.pred.test, as.numeric(y)), 1, function(x) {mean(x[-(Bsample + 1)] != x[Bsample + 1],na.rm = TRUE)}), na.rm = TRUE)
	#cat("\n  one :",one,"\n") 

	resubst=mean(err.inbag, na.rm=T)
	#resubst <- mean(mat.pred.inbag != as.numeric(y), na.rm=T)
	#cat("\n  resubst :",resubst,"\n") 

	err632 <- 0.368 * resubst + 0.632 * one

	gamma <-sum(outer(as.numeric(y),as.numeric(mat.pred.inbag),function(x, y) ifelse(x == y, 0, 1)),na.rm=T)/(length(y)^2)
	
	r <- (one - resubst)/(gamma - resubst)
	r <- ifelse(one > resubst & gamma > resubst, r, 0)
	if((r > 1) | (r < 0)) { ## just debugging; 
        #print(paste("r outside of 0, 1 bounds: one", one,
         #           "resubst", resubst, "gamma", gamma))
        if(r > 1) {
            r <- 1
            #print("setting r to 1")
        }
        else if(r < 0) {
            r <- 0
            #print("setting r to 0")
        }
    	}
    
    	errprime <- min(one, gamma)
    	err <- err632 + (errprime - resubst) *(0.368 * 0.632 * r)/(1 - 0.368 * r)

	##fill the error vector
	error.boot[ngene]=err

    } #end loop on ngene


	out <-list(
	   maxvar=maxvar,
           nvar=nvar,
	   weight.eval=weight,
           ntreeTest=ntreeTest,	
  	   matTrain=obj$matTrain,
	   matProb=obj$matProb,
	   error=error.boot,
	   sampleWeight= if(!weight) NULL else {sampleWeight},
	   matPredInbag=mat.pred.inbag,
	   matPredTest=mat.pred.test
	)

	return(out)


  }
