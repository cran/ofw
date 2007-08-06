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

"evaluateSVM.learnSVM"<-
  function(obj,
           maxvar=15,
           weight=F,
	   ...) {
    
    x=obj$x
    y=obj$y
    nsample=obj$nsample
    nvariable=obj$nvariable
    niteration=obj$niteration
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

    if(obj$weight.learn==T  && weight==T){
    sampleWeight=obj$weightingOption$sampleWeight  #if obj has been weighted
    classWeight=obj$weightingOption$classWeight
    }	

    #check input
    if(maxvar >=nvariable) stop("maxvar shoud not be greater than the number of variables")

    ##compute the weights (if needed) for each sample
    if(weight==T && obj$weight.learn==F){
	for(boot in 1:Bsample){
		numWeight[,boot]=summary(y[obj$matTrain[,boot]])
		classWeight[,boot]=1/(numWeight[,boot] * nlevels(y[obj$matTrain[,boot]]))    
		#classWeight[,boot]=nsample/numWeight[,boot]
		#classWeight[,boot]=classWeight[,boot]/sum(classWeight[,boot])
		for(n in 1:nlevels(y)){
		sampleWeight[which(as.integer(y[obj$matTrain[,boot]])==n),boot]= classWeight[n,boot]
	}
	}

    }
	
############################### EVAL ERREUR ############################################
eval.erreur=function(x,y,train,test, P,ngene){
	liste=as.numeric(names(sort(P, decreasing=T)[1:ngene]))
	data.train=x[train,liste]     
	data.test=x[test,liste]

	##pour l'erreur de resubstitution
	svm.train=svm(data.train, y[train],kernel="linear")
	sample.pred.train[train]=predict(svm.train, data.train)
	
	##pour erreur interne oob
	sample.pred.test[test]=predict(svm.train, data.test)
	
	return(list(sample.pred.train[train],sample.pred.test[test]))
} #fin eval.erreur
###############################   EVAL ERREUR WEIGHT    ################################
#evaluation de l erreur pour calculer e632+
eval.erreur.weight=function(x,y,train,test, P,ngene, vect.w.train, vect.w.test){
	liste=as.numeric(names(sort(P, decreasing=T)[1:ngene]))
	data.train=x[train,liste]     
	data.test=x[test,liste]


	##pour l'erreur de resubstitution
	svm.train=svm(data.train, y[train],kernel="linear")
	sample.pred.train[train]=predict(svm.train, data.train)
	pred.train=predict(svm.train, data.train)
	err.train = sum(vect.w.train[which(pred.train != y[train])])/sum(vect.w.train)  #c bien la somme?


	##pour erreur interne oob
	sample.pred.test[test]=predict(svm.train, data.test)
	pred.test=predict(svm.train, data.test)
	err.test = sum(vect.w.test[which(pred.test != y[test])])/sum(vect.w.test)



	return(list(sample.pred.train[train],sample.pred.test[test], err.train, err.test))
} #fin eval.erreur.weight

###########################################################################
#pr l erreur .632+
#erreur de resubstitution
sample.pred.train=vector(length=nsample)

#erreur in sample
sample.pred.test=vector(length=nsample)

error.boot=vector(length=maxvar)              #erreur e632+

mat.pred.inbag=matrix(nrow=nsample, ncol=Bsample)
mat.pred.test=matrix(nrow=nsample, ncol=Bsample)
if (weight==T) {
	err.inbag=vector(length=Bsample, mode="numeric")
	err.test=vector(length=Bsample, mode="numeric")
}

for (ngene in 1:maxvar)
{
#print(ngene)

for (boot in 1:Bsample)
{
#pr chaque echantillon bootstrap on travaille sur training et test set boot

train.b = obj$matTrain[,boot]
test.b = setdiff(1:nsample, train.b)
P.b=obj$matProb[,boot]
names(P.b)=c(1:nvariable)
if (weight==T) {
	vect.w.train=sampleWeight[train.b,boot]
	vect.w.test=sampleWeight[test.b,boot]
}


### Call the function and collect the results
if (weight==T) {
	res.erreur=eval.erreur.weight(x=x, y=y,train=train.b, test=test.b, P=P.b,ngene=ngene,vect.w.train=vect.w.train, vect.w.test=vect.w.test)
	} else {
	res.erreur=eval.erreur(x=x, y=y,train=train.b, test=test.b, P=P.b,ngene=ngene)
}

mat.pred.inbag[train.b,boot]=res.erreur[[1]]
mat.pred.test[test.b,boot]=res.erreur[[2]]
if (weight==T) {
err.inbag[boot]=res.erreur[[3]]
err.test[boot]=res.erreur[[4]]
}

} #fin boucle boot

	#ce code e632+ provient du code en R thorsten ds le package ipred ainsi que de Diaz ds package
	#varSelRF. J ai verifie les fonctions c est ok

# c est l erreur leave one out pr chaque ech bootstrap
if (weight==T) { 
	one=mean(err.test, na.rm=T) 
	} else { 
	one <- mean(apply(cbind(mat.pred.test, as.numeric(y)), 1, function(x) {mean(x[-(Bsample + 1)] != x[Bsample + 1],na.rm = TRUE)}), na.rm = TRUE) 
	}

if(weight==T) {
	resubst=mean(err.inbag, na.rm=T)
	} else {
	resubst <- mean(mat.pred.inbag != as.numeric(y), na.rm=T)
	}

err632 <- 0.368 * resubst + 0.632 * one
gamma <-sum(outer(as.numeric(y),as.numeric(mat.pred.inbag),function(x, y) ifelse(x == y, 0, 1)),na.rm=T)/(length(y)^2)


r <- (one - resubst)/(gamma - resubst)
    r <- ifelse(one > resubst & gamma > resubst, r, 0)
    if((r > 1) | (r < 0)) { ## just debugging; eliminar mÃ¡s adelante
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


#print(error.boot[ngene])

} #fin boucle erreur sur tous les ngene

##############################################################

	out <-list(
	   maxvar=maxvar,
	   weight.eval=weight,
	   #niteration=niteration,
  	   matTrain=obj$matTrain,
	   matProb=obj$matProb,
	   error=error.boot,
	   sampleWeight=sampleWeight,
	   matPredInbag=mat.pred.inbag,
	   matPredTest=mat.pred.test
	)

	return(out)


  }
