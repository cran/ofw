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

"learnCART.default"<-
   function(x, 
	    y,
	    ntree=50, 
	    nforest=10, 
	    mtry=5, 
	    do.trace=FALSE,
            nstable=50,
            weight=FALSE,
            Bsample=5, ...) {


	mat.train=matrix(nrow=nrow(x), ncol=Bsample)
	mat.prob=matrix(nrow=ncol(x), ncol=Bsample)
	mat.pred.inbag=matrix(nrow=nrow(x), ncol=Bsample)
	mat.pred.test=matrix(nrow=nrow(x), ncol=Bsample)

	#for the weights if needed
	if(weight==T){
	numWeight=matrix(nrow=nlevels(y), ncol=Bsample)
	classWeight=matrix(nrow=nlevels(y),ncol=Bsample)
	sampleWeight=matrix(nrow=nrow(x),ncol=Bsample)
	}
	
	nsample=nrow(x)
	nvariable=ncol(x)
	
	
	if (!is.factor(y)) stop("y is not a factor")
	if (length(unique(y)) < 2) stop("Need at least two classes to do classification.")

	## Make sure mtry is in reasonable range.
	if (mtry < 1 || mtry > nvariable) warning("invalid mtry: reset to within valid range")
	#mtry <- max(1, min(p, round(mtry)))

	if (length(y) != nsample) stop("length of response must be the same as predictors")

	## Check for NAs.
	if (any(is.na(x))) stop("NA not permitted in predictors")
	if (any(is.na(y))) stop("NA not permitted in response")

	## Check for empty classes:
	if (any(table(y) == 0)) stop("Can't have empty classes in y.")



	#learn the Bsample algorithms on the Bsample samples
	for (boot in 1:Bsample)
	{
	##in case one sample is empty in the bootstrap:
	cont=T
	while (cont==T){
	train=sample(1:nsample, nsample, replace=T)
	if (any(table(y[train]) < min(table(y)))) {cont=T} else {cont=F}
	}
 
	mat.train[,boot]=train
	cat("\n  Learning ofwCART on the boostrap sample",boot,"\n")
	xprim=x[train,]
	yprim=y[train]
	obj=ofwCART(xprim,yprim, ntree=ntree, nforest=nforest, mtry=mtry, weight=weight, do.trace=do.trace, nstable=nstable)
	mat.prob[,boot]=obj$prob

	#compute the weights for each sample
	if(weight==T){
		numWeight[,boot]=summary(y[train])
		classWeight[,boot]=nsample/numWeight[,boot]
		classWeight[,boot]=classWeight[,boot]/sum(classWeight[,boot])
		for(n in 1:nlevels(y)){
		sampleWeight[which(as.integer(y[train])==n),boot]= classWeight[n,boot]
		}
	} #fin weight
	}#fin boot




	out <-list(x=x,
                   y=y,
                   nsample=nsample,
	           nclass=nlevels(y),
                   nvariable=nvariable,
                   weight.learn=weight,
		   Bsample=Bsample,
  		   matTrain=mat.train,
		   matProb=mat.prob,
		   weightingOption= if (!weight) NULL else {list(classWeight=classWeight, sampleWeight=sampleWeight) })

        class(out) <- "learnCART"

	out
}
