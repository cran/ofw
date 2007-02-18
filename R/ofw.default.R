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

"ofw" <- function(x, ...) UseMethod("ofw")

"ofw.default"<-
   function(    x,
   	    	y,
   	    	type="CART",
   	    	ntree= if(type=="CART") 50 else NULL,
   	    	nforest= if(type=="CART") 100 else NULL,
   	    	nsvm= if(type=="SVM") 100 else NULL,
   	    	mtry=5,
   	    	do.trace=FALSE,
   	    	nstable=25,
   	    	keep.inbag=if(type=="CART")  FALSE else NULL,
   	    	keep.forest=if(type=="CART") TRUE else NULL,
   	    	weight=FALSE, ...){


 
#pour cart:
#function(x, y, ntree=50, nforest=100, mtry=5, do.trace=FALSE, nstable=25, keep.inbag=FALSE, keep.forest=TRUE, weight=FALSE)

if (type=="CART"){

	ofw.cart=ofwCART(x, y, ntree=ntree, nforest=nforest, mtry=mtry, do.trace=do.trace, nstable=nstable, keep.inbag=keep.inbag, keep.forest=keep.forest, weight=weight)

	cl <- match.call()
        cl[[1]] <- as.name("ofw")
        out <- list(call = cl,
        	    type=type,
                    classes = levels(y),
	            mean.error=ofw.cart$mean.error,
		    prob=ofw.cart$prob,      
		    ntree = ntree,
		    nforest=nforest,
		    maxiter=ofw.cart$maxiter, 
                    mtry = mtry,
	            do.trace=do.trace,
		    nstable=nstable,
		    weight=weight,
		    weightingOption= if (!weight) NULL else {list(classWeight=ofw.cart$weightingOption$classWeight, sampleWeight=ofw.cart$weightingOption$sampleWeight) },
                    forest = ofw.cart$forest,
                    inbag = ofw.cart$inbag)
}
            
#pour svm:
#function(x, y, nsvm=100, mtry=5, do.trace=FALSE, nstable=25, weight=FALSE)
if (type=="SVM"){

	ofw.svm=ofwSVM(x, y, nsvm=nsvm, mtry=mtry, do.trace=do.trace, nstable=nstable, weight=weight)
	
	cl <- match.call()
        cl[[1]] <- as.name("ofw")

	out <- list(
		call = cl,
		type=type,
		classes = levels(y),
	    	prob=ofw.svm$prob,
		nsvm=nsvm,
		maxiter=ofw.svm$maxiter, 
		mtry = mtry,
		do.trace=do.trace,
		nstable=nstable,
		weight=weight,
		weightingOption= if (!weight) NULL else {list(classWeight=ofw.svm$weightingOption$classWeight, sampleWeight=ofw.svm$weightingOption$sampleWeight) }
                )
}
	class(out) <- "ofw"
        return(out)
        

}
