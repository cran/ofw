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


"ofwSVM" <-function(x, ...)  UseMethod("ofwSVM")

"ofwSVM.default" <-
    function(
	     x, 
	     y, 
	     nsvm=100,
             mtry=5, 
	     do.trace=FALSE,
	     nstable=25,
	     weight=FALSE,
	      ...) {

nsample=nrow(x)
nvariable=ncol(x)

if (!is.factor(y)) stop("y is not a factor")
if (length(unique(y)) < 2) stop("Need at least two classes to do classification.")

## Make sure mtry is in reasonable range.
if (mtry < 1 || mtry > nvariable) warning("invalid mtry: reset to within valid range")

if (length(y) != nsample) stop("length of response must be the same as predictors")

## Check for NAs.
if (any(is.na(x))) stop("NA not permitted in predictors")
if (any(is.na(y))) stop("NA not permitted in response")

## Check for empty classes or small classes:
if (any(table(y) == 0)) stop("Can't have empty classes in y.")
#if (any(table(y) == 2)) stop("Can't have less than 3 cases per class in y.")

x = as.data.frame(x)


##########################  MAIN   ####################################
#compute the weights for each sample
if(weight){
	#declarations
	numWeight=vector(length=nlevels(y))
	classWeight=vector(length=nlevels(y))
	sampleWeight=vector(length=nsample)

	numWeight=summary(y)
	classWeight=1/(numWeight * nlevels(y))                    
	for(n in 1:nlevels(y)){
	sampleWeight[which(as.integer(y)==n)]= classWeight[n]
	}
} #fin weight

#declarations
P = c(rep(1/nvariable, nvariable))	 #initialisation : P est un vecteur uniforme
names(P)=colnames(x)
W= vector(length=mtry)                     #sous ensemble W de variables de longeur mtry
erreur.it= vector(mode="numeric", length=nsvm)  #erreur moyenne de chaque iteration
G.final=vector(length=length(P), mode="numeric")   #gradient de l energie

j=1
iter = 1
iter.stop=0

while(iter <=nsvm)
{
#learn the svm on a bootstrap sample
res.svm<-svmlearn(x=x, y=y,P=P, lw=mtry, nsample = nsample, nvariable = nvariable, vectWeight = classWeight, weight = weight)

#get the gradient G
G.final=res.svm[[1]]

#and the OOB error rate
erreur.it[iter]=res.svm[[2]]

P=P - G.final/(1+iter)               #gradient descent  
P=simplexe(P)                        #projection on the simplex

# in case of early stopping criterion
if(do.trace){
if (iter==(j*do.trace))
	{
		if (j==1) {l1=names(sort(P, decreasing=TRUE)[1:nstable])} else{l2=names(sort(P, decreasing=TRUE)[1:nstable])}
		cat("\n", "iteration", iter, "     " )
		if (j >1) cat("stable variables: ", length(intersect(l1, l2)), " ")
		if (j>1) {if (length(intersect(l1, l2))==nstable) {iter.stop=iter;iter=nsvm;} else {l1=l2}}
		j=j+1
	} 
}

iter=iter+1
}  #end for nsvm
cat("\n")

if(iter.stop != 0) {maxiter=iter.stop} else {maxiter=nsvm}

	cl <- match.call()
        cl[[1]] <- as.name("ofwSVM")
out <- list(
		call = cl,
		type = "classification", 
		classes = levels(y),
	    	prob=P,
		nsvm=nsvm,
		maxiter=maxiter, 
		mtry = mtry,
		do.trace=do.trace,
		nstable=nstable,
		weight=weight,
		classWeight= if (!weight) NULL else classWeight, 
		sampleWeight=  if (!weight) NULL else sampleWeight
                )
                   

    class(out) <- "ofwSVM"
    return(out)



} #fin ofwSVM.default
