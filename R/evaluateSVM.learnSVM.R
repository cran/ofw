#
# Copyright (C) 2007 Kim-Anh L� Cao, Patrick Chabrier, INRA,
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

"evaluateSVM" <- function(obj, ...) UseMethod("evaluateSVM")

"evaluateSVM.learnSVM"<-
  function(obj,
           maxvar=15,
           weight=FALSE,
	   ...) {
    
    x=obj$x
    y=obj$y
    nsample=obj$nsample
    nvariable=obj$nvariable
    nsvm=obj$nsvm
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

    if(obj$weight && weight){
    sampleWeight=obj$sampleWeight  #if obj has been weighted
    classWeight=obj$classWeight
    }	

    #check input
    if(maxvar >=nvariable) stop("maxvar shoud not be greater than the number of variables")

    ##compute the weights (if needed) for each sample
    if(weight && !obj$weight){
	for(boot in 1:Bsample){
		numWeight[,boot]=summary(y[obj$matTrain[,boot]])
		classWeight[,boot]=1/(numWeight[,boot] * nlevels(y[obj$matTrain[,boot]]))    
		for(n in 1:nlevels(y)){
		sampleWeight[which(as.integer(y[obj$matTrain[,boot]])==n),boot]= classWeight[n,boot]
	}
	}

    }
	
###############################   EVAL ERREUR    ################################
#evaluation de l erreur pour calculer e632+
eval.erreur=function(x,y,train,test, P,ngene, vect.w.train=NULL, vect.w.test=NULL, weight){
	liste=as.numeric(names(sort(P, decreasing=TRUE)[1:ngene]))
	data.train=x[train,liste]     
	data.test=x[test,liste]


	##pour l'erreur de resubstitution
	svm.train=svm(data.train, y[train],kernel="linear")
	sample.pred.train[train]=predict(svm.train, data.train)
	if (weight) {
	err.train = sum(vect.w.train[which(sample.pred.train[train] != y[train])])/sum(vect.w.train)
	} else {
	err.train = length(which(sample.pred.train[train] != y[train]))/length(train)
	}

	##pour erreur interne oob
	sample.pred.test[test]=predict(svm.train, data.test)
	if (weight) {
	err.test = sum(vect.w.test[which(sample.pred.test[test] != y[test])])/sum(vect.w.test)
	} else {
	err.test = length(which(sample.pred.test[test] != y[test]))/length(test)
	}

	return(list(sample.pred.train[train],sample.pred.test[test], err.train, err.test))
} #fin eval.erreur

###########################################################################
#pr l erreur .632+
#erreur de resubstitution
sample.pred.train=vector(length=nsample)

#erreur in sample
sample.pred.test=vector(length=nsample)

error.boot=vector(length=maxvar)              #erreur e632+

mat.pred.inbag=matrix(nrow=nsample, ncol=Bsample)
mat.pred.test=matrix(nrow=nsample, ncol=Bsample)
err.inbag=vector(length=Bsample, mode="numeric")
err.test=vector(length=Bsample, mode="numeric")


    ##BIG LOOP ON MAXvar
cat("\n  Calculating e632+ for each of the ",maxvar," variables \n")
for(ngene in 1:maxvar){
	cat("\n  variable :",ngene,"\n") 

    ##LOOP on the bootstrap samples
	for (boot in 1:Bsample)
	{
	#pr chaque echantillon bootstrap on travaille sur training et test set boot

	train.b = obj$matTrain[,boot]
	test.b = setdiff(1:nsample, train.b)
	P.b=obj$matProb[,boot]
	names(P.b)=c(1:nvariable)

	### Call the function and collect the results
	if (weight) {
	res.erreur=eval.erreur(x=x, y=y,train=train.b, test=test.b, P=P.b,ngene=ngene,vect.w.train=sampleWeight[train.b,boot], vect.w.test=sampleWeight[test.b,boot], weight=weight)
	} else {
	res.erreur=eval.erreur(x=x, y=y,train=train.b, test=test.b, P=P.b,ngene=ngene, weight=weight)
	}

mat.pred.inbag[train.b,boot]=res.erreur[[1]]
mat.pred.test[test.b,boot]=res.erreur[[2]]
err.inbag[boot]=res.erreur[[3]]
err.test[boot]=res.erreur[[4]]

} #fin boucle boot

	#ce code e632+ provient du code en R thorsten ds le package ipred ainsi que de Diaz ds package
	#varSelRF. J ai verifie les fonctions c est ok

	# erreur leave one out pr chaque ech bootstrap
	one=mean(err.test, na.rm=TRUE) 

	resubst=mean(err.inbag, na.rm=TRUE)


	err632 <- 0.368 * resubst + 0.632 * one
	gamma <-sum(outer(as.numeric(y),as.numeric(mat.pred.inbag),function(x, y) ifelse(x == y, 0, 1)),na.rm=TRUE)/(length(y)^2)

	r <- (one - resubst)/(gamma - resubst)
    	r <- ifelse(one > resubst & gamma > resubst, r, 0)
    	if((r > 1) | (r < 0)) { ## just debugging;
        print(paste("r outside of 0, 1 bounds: one", one,
                    "resubst", resubst, "gamma", gamma))
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

	error.boot[ngene]=err



} #fin boucle erreur sur tous les ngene

##############################################################

	out <-list(
	   maxvar=maxvar,
	   weight.eval=weight,
  	   matTrain=obj$matTrain,
	   matProb=obj$matProb,
	   error=error.boot,
	   sampleWeight=sampleWeight,
	   matPredInbag=mat.pred.inbag,
	   matPredTest=mat.pred.test
	)

	return(out)


  }
