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

"ofwTune" <- function(x, ...) UseMethod("ofwTune")
  
  
"ofwTune.default" <-
    function(
	     x, 
	     y,
	     type="CART",
	     ntree= if(type=="CART") 50 else NULL, 
	     nforest= if(type=="CART") 100 else NULL,
	     nsvm= if(type=="SVM") 500 else NULL,
             mtry.test=seq(5,15,length=3), 
	     do.trace=FALSE,
	     nstable=10,
             weight=FALSE, ...) {
             
             
	param=matrix(nrow=2, ncol=length(mtry.test))
	rownames(param)=c("mtry", "length")
	colnames(param)=c(1:length(mtry.test))
	iter.max=matrix(nrow=2 , ncol=length(mtry.test))
	colnames(iter.max)=c(1:length(mtry.test))
	

	if (type=="CART"){
		rownames(iter.max)=c("ofwCART1", "ofwCART2")
		for(i in 1:length(mtry.test)){
		mtry=mtry.test[i]
		cat("Computing ofwCART1 and ofwCART2 for mtry=", mtry, ": \n",sep="")
		cat("ofwCART1", "\n")
		ofwCART1=ofwCART(x, y, ntree=ntree, nforest=nforest, mtry=mtry, nstable=nstable, do.trace=do.trace, weight=weight)
		iter.max[1,i]=ofwCART1$maxiter
		cat("ofwCART2","\n")
		ofwCART2=ofwCART(x, y,  ntree=ntree, nforest=nforest, mtry=mtry, nstable=nstable, do.trace=do.trace, weight=weight)
		iter.max[2,i]=ofwCART2$maxiter
		l=length(intersect(names(sort(ofwCART1$prob, decreasing=TRUE)[1:nstable]),names(sort(ofwCART2$prob, decreasing=TRUE)[1:nstable])))
		param[,i]=c(mtry, l)
		}
	}
	
	if (type=="SVM"){
		rownames(iter.max)=c("ofwSVM1", "ofwSVM2")
		for(i in 1:length(mtry.test)){
		mtry=mtry.test[i]
		cat("Computing ofwSVM1 and ofwSVM2 for mtry=", mtry, ": \n",sep="")
		cat("ofwSVM1", "\n")
		ofwSVM1=ofwSVM(x, y, nsvm=nsvm, mtry=mtry, nstable=nstable, do.trace=do.trace, weight=weight)
		iter.max[1,i]=ofwSVM1$maxiter
		cat("ofwSVM2", "\n")
		ofwSVM2=ofwSVM(x, y, nsvm=nsvm, mtry=mtry, nstable=nstable, do.trace=do.trace, weight=weight)
		iter.max[2,i]=ofwSVM2$maxiter
		l=length(intersect(names(sort(ofwSVM1$prob, decreasing=TRUE)[1:nstable]),names(sort(ofwSVM2$prob, decreasing=TRUE)[1:nstable])))
		param[,i]=c(mtry, l)
		}
	}
		
		
	cl <- match.call()
        cl[[1]] <- as.name("ofwTune")	
			

	out=list(
		call = cl,
		type=type,
		nstable=nstable,
		do.trace=do.trace,
		mtry.test=mtry.test,
		param=param,
		itermax=iter.max,
		iteration= if (is.null(nforest)) nsvm else nforest,
		weight=weight)
	class(out) <- "ofwTune"
	return(out)

}



