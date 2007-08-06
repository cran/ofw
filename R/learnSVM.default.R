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

"learnSVM.default" <-
   function(x, 
	    y, 
	    niteration=1000, 
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

## Check for empty classes:
if (any(table(y) == 0)) stop("Can't have empty classes in y.")




################FONCTIONS ######
simplexe=function(P)
{
#1 Projeter sur l'hyperplan Hf si P n'appartient pas a Hf
somme=sum(P)
lambda=(somme-1)/length(P)
P=P-lambda  #projection sur Hf
#2Projeter sur le simplexe les elements de P <0
Sf= 1		#pas sur le simplexe (valeur arbitraire !=0)
while(Sf!=0)
{
Sf=0
J=sum(P<=0)          #nb d'elements de P <=0
somme=sum(P[P>0])    #somme des elements de P >0
P[P<=0]=0            #on met les proba negatives a zero
P[P>0]=P[P>0] + (1-somme)/(length(P) - J)   #projection surle simplexe
Sf=sum(P<0)    #y a t il encore des elements hors du simplexe ?
}
return(P)
}
######################## SVMSLAVE ###################################
svmslave=function(xprim, yprim,lw, P,trainsample) 
{
G=vector(length=length(P))
W=sample(nvariable, size=lw, prob=P, replace=F) #tirer W selon la proba P

cont=T
while (cont==T){
train=sample(trainsample, nsample, replace=T)
if ((any(table(yprim[train]) == 0))|| (any(table(yprim[setdiff(trainsample, train)]) == 0))) {cont=T} else {cont=F}
}
test = setdiff(trainsample, train)                 #echantillon oob     
data.train = xprim[train,W]       
data.test=xprim[test,W]
svm.train=svm(data.train, yprim[train], kernel="linear")
mat=table(yprim[test],predict(svm.train, data.test))
erreur= (length(test) - sum(diag(mat)))/length(test)
G[W]= G[W] +erreur/(1000*P[W])  
return(list(G,erreur))
}
###########################  SVMSLAVEWEIGHT  ######################################
svmslaveWeight=function(xprim, yprim,lw, P, vectWeight, trainsample) 
{
G=vector(length=length(P), mode="numeric")
W=sample(nvariable, size=lw, prob=P, replace=F) #tirer W selon la proba P
cont=T
while (cont==T){
train=sample(trainsample, nsample, replace=T)
if ((any(table(yprim[train]) == 0))|| (any(table(yprim[setdiff(trainsample, train)]) == 0))) {cont=T} else {cont=F}
}
test = setdiff(trainsample, train)                 #echantillon oob     
data.train = xprim[train,W]       
data.test=xprim[test,W]
svm.train=svm(data.train, yprim[train], kernel="linear")
mat=table(yprim[test],predict(svm.train, data.test))
erreur=sum((apply(mat, 1, sum) -diag(mat))*classWeight)/sum(apply(mat, 1, sum)*classWeight)
G[W]= G[W] +erreur/(1000*P[W])  
return(list(G,erreur))
}



##########################  MAIN   ####################################

#for the weights if needed
if(weight==T){  
numWeight=matrix(nrow=nlevels(y), ncol=Bsample)
classWeight=matrix(nrow=nlevels(y), ncol=Bsample)
sampleWeight=matrix(nrow=nsample,ncol=Bsample)
}

mat.train=matrix(nrow=nsample, ncol=Bsample)
matrice.P=matrix(nrow=nvariable, ncol=Bsample)




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

	#compute the weights for each sample
	if(weight==T){
		numWeight[,boot]=summary(y[train])
		classWeight[,boot]=1/(numWeight[,boot] * nlevels(y[train]))    
		#classWeight[,boot]=nsample/numWeight[,boot]
		#classWeight[,boot]=classWeight[,boot]/sum(classWeight[,boot])
		for(n in 1:nlevels(y)){
		sampleWeight[which(as.integer(y[train])==n),boot]= classWeight[n,boot]
		}
	} #fin weight
	}#fin boot




for (boot in 1:Bsample)
{
cat("\n  Learning ofwSVM on the boostrap sample",boot,"\n")

train=mat.train[,boot]
xprim=x[train,]
yprim=y[train]


P = c(rep(1/nvariable, nvariable))	 #initialisation : P est un vecteur uniforme
W= vector(length=mtry)  #sous ensemble W de variables de longeur mtry
erreur.it= vector(mode="numeric", length=niteration)  #erreur moyenne de chaque iteration
G.final=vector(length=length(P), mode="numeric")   #gradient de l energie

j=1
iter = 1
iter.stop=0

while(iter <= niteration)
{
# Call the function in all the children, and collect the results
if(weight==T) {res.slaves <-svmslaveWeight(xprim, yprim,P=P, lw=mtry, vectWeight=classWeight[,boot], trainsample=mat.train[,boot])}  else {res.slaves <-svmslave(xprim, yprim,P=P, lw=mtry, trainsample=mat.train[,boot])}

#on recupere le gradient G et l'erreur
G.final=res.slaves[[1]]
erreur.it[iter]=res.slaves[[2]]


P=P - G.final/(10+iter)               #descente de gradient  
P=simplexe(P)                        #projection sur le simplexe 


if(do.trace){
if (iter==(j*do.trace))
{
if (j==1) {l1=which(P>=(sort(P, decreasing=T)[nstable]))} else{l2=which(P>=(sort(P, decreasing=T)[nstable]))}
if (j>1) {if (length(intersect(l1, l2))==nstable) {iter.stop=iter;iter=niteration;} else {l1=l2}}
j=j+1
} 
}

iter=iter+1

}  #fin boucle niteration#if(iter.stop != 0) {cat("iteration max",iter.stop,"\n")} else {cat("iteration max",niteration,"\n")}

matrice.P[,boot]=P

}  #fin boot
####a ce stade on a appris les lois P_{boot}

########################################################################"
	out <-list(x=x,
                   y=y,
                   nsample=nsample,
	           nclass=nlevels(y),
                   nvariable=nvariable,
                   weight.learn=weight,
		   niteration=niteration,
		   Bsample=Bsample,
		   #mtry=mtry,
  		   matTrain=mat.train,
		   matProb=matrice.P, 
		   weightingOption= if (!weight) NULL else {list(classWeight=classWeight, sampleWeight=sampleWeight) }
		   )

        class(out) <- "learnSVM"
	return(out)
} #fin learnSVM.default


