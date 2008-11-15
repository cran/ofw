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

"learnSVM" <- function(x, ...) UseMethod("learnSVM")

"learnSVM.default" <-
   function(x, 
	    y, 
	    nsvm=1000, 
	    mtry=10, 
	    do.trace=FALSE,
            nstable=50,
            weight=FALSE,
            Bsample=5,...) {


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

## Check for empty classes or small classes:
if (any(table(y) == 0)) stop("Can't have empty classes in y.")
if (any(table(y) == 2)) stop("Can't have less than 3 cases per class in y for evaluation.")


##########################  MAIN   ####################################

#for the weights if needed
if(weight){  
	classWeight=matrix(nrow=nlevels(y), ncol=Bsample)
	sampleWeight=matrix(nrow=nsample,ncol=Bsample)
	}

mat.train=matrix(nrow=nsample, ncol=Bsample)
matrice.P=matrix(nrow=nvariable, ncol=Bsample)


#learn the Bsample algorithms on the Bsample samples
	for (boot in 1:Bsample)
	{
	#in case one sample is empty in the bootstrap:
	cont=TRUE
	while (cont){
		train=sample(1:nsample, nsample, replace=TRUE)
		if (any(table(y[train]) < 2)) {cont=TRUE} else {cont=FALSE}
		}
	mat.train[,boot]=train
	}#fin boot

# loop on the bootstrap samples
for (boot in 1:Bsample)
{
cat("\n", " Learning ofwSVM on the boostrap sample",boot,"\n")

train=mat.train[,boot]
xprim=x[train,]
yprim=y[train]

boot.svm = ofwSVM( x = xprim, y = yprim, nsvm=nsvm, mtry=mtry, do.trace=do.trace, nstable = nstable, weight=weight)

matrice.P[,boot]=boot.svm$prob

if (weight) {
	classWeight[,boot] = boot.svm$classWeight
	sampleWeight[,boot]= boot.svm$sampleWeight
	}

}  #fin boot
cat("\n")
#a ce stade on a appris les lois P_{boot}

########################################################################"
	out <-list(x=x,
                   y=y,
                   nsample=nsample,
	           nclass=nlevels(y),
                   nvariable=nvariable,
                   nsvm=nsvm,
		   Bsample=Bsample,
  		   matTrain=mat.train,
		   matProb=matrice.P, 
		   weight=weight,
		   classWeight= if (!weight) NULL else classWeight, 
		   sampleWeight=  if (!weight) NULL else sampleWeight
		   )

        class(out) <- "learnSVM"
	return(out)
} #fin learnSVM.default


