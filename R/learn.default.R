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

"learn" <- function(x, ...) UseMethod("learn")

"learn.default"<-
   function(x, 
	    y,
	    type="CART",
	    ntree= if(type=="CART") 50 else NULL, 
	    nforest= if(type=="CART") 100 else NULL,
	    nsvm= if(type=="SVM") 20000 else NULL,
            mtry=5, 
	    do.trace=FALSE,
            nstable=50,
            weight=FALSE,
            Bsample=5, ...) {
#pour cart:
#function(x, y, ntree=50,  nforest=10, mtry=5, do.trace=FALSE, nstable=50, weight=FALSE, Bsample=5)
if (type=="CART"){

#faire des warnings

	learn.cart=learnCART(x, y, ntree=ntree,  nforest=nforest, mtry=mtry, do.trace=do.trace, nstable=nstable, weight=weight, Bsample=Bsample)


	out <-list(x=x,
                   y=y,
                   type=type,
                   nsample=learn.cart$nsample,
	           nclass=nlevels(y),
                   nvariable=learn.cart$nvariable,
                   Bsample=Bsample,
  		   matTrain=learn.cart$matTrain,
		   matProb=learn.cart$matProb,
		   weight=weight,
		   classWeight= if (!weight) NULL else learn.cart$classWeight, 
		   sampleWeight=  if (!weight) NULL else learn.cart$sampleWeight 
		   )
		   }
            
#pour svm:
#function(x, y, nsvm=1000, mtry=10, do.trace=FALSE, nstable=50, weight=FALSE, Bsample=5)
if (type=="SVM"){

	learn.svm=learnSVM(x, y, nsvm=nsvm, mtry=mtry, do.trace=do.trace, nstable=nstable, weight=weight, Bsample=Bsample)

	out <-list(x=x,
                   y=y,
                   type=type,
                   nsample=learn.svm$nsample,
	           nclass=nlevels(y),
                   nvariable=learn.svm$nvariable,
                   nsvm=learn.svm$nsvm,
		   Bsample=Bsample,
  		   matTrain=learn.svm$matTrain,
		   matProb=learn.svm$matProb,
		   weight=weight, 
		   classWeight= if (!weight) NULL else learn.svm$classWeight, 
		   sampleWeight=  if (!weight) NULL else learn.svm$sampleWeight 
		   )
		   
}

        class(out) <- "learn"
        return(out)

}
