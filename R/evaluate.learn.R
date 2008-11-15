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

"evaluate" <- function(obj, ...) UseMethod("evaluate")


"evaluate.learn"<-
   function(obj, 
	    maxvar=15,
	    type=obj$type,
	    nvar=if(obj$type=="CART")  obj$nclass+1 else NULL,
	    ntreeTest=if(obj$type=="CART") 100 else NULL,
	    weight=FALSE,...) {

 
#pour cart:
#function(obj, maxvar=15, nvar=obj$nclass+1, ntreeTest=100, weight=FALSE)
if (type=="CART"){

#faire des warnings

	eval.cart=evaluateCART.learnCART(obj, maxvar=maxvar, nvar=nvar, ntreeTest=ntreeTest, weight=weight)


	out <-list(
	   maxvar=maxvar,
           nvar=nvar,
	   weight.eval=weight,
	   weight.learn=obj$weight.learn,
           ntreeTest=ntreeTest,	
  	   matTrain=obj$matTrain,
	   matProb=obj$matProb,
	   error=eval.cart$error,
	   sampleWeight= if(!weight) NULL else {eval.cart$sampleWeight},
	   matPredInbag=eval.cart$matPredInbag,
	   matPredTest=eval.cart$matPredTest
	)
}
            
if (type=="SVM"){

	eval.svm=evaluateSVM.learnSVM(obj, maxvar=maxvar, nvar=nvar, weight=weight)

	out <-list(
	   maxvar=maxvar,
	   weight.eval=weight,
	   weight.learn=obj$weight.learn,
  	   matTrain=obj$matTrain,
	   matProb=obj$matProb,
	   error=eval.svm$error,
	   sampleWeight=eval.svm$sampleWeight,
	   matPredInbag=eval.svm$matPredInbag,
	   matPredTest=eval.svm$matPredTest
	)
}

        return(out)

}
