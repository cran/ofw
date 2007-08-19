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

"ofwCART" <- function(x, ...) UseMethod("ofwCART")
  
  
"ofwCART.default" <-
    function(
	     x, 
	     y, 
	     ntree=50, 
	     nforest=100,
             mtry=5, 
	     do.trace=FALSE,
	     nstable=25,
             keep.inbag=FALSE,
	     keep.forest=TRUE,
	     weight=FALSE, ...) {

    classRF <- is.factor(y)
	if (!classRF){warning("Regression is not possible, check that y is a factor")}	

	if (classRF  && length(unique(y)) < 2)
        stop("Need at least two classes to do classification.")
        
        ##if (do.trace==TRUE) do.trace=nforest
    if (!is.factor(y))
	stop("y is not a factor!")
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    keep.forest <- keep.forest
    
    ## Make sure mtry is in reasonable range.
    if (mtry < 1 || mtry > p) warning("invalid mtry: reset to within valid range")

    if (length(y) != n) stop("length of response must be the same as predictors")


    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (any(is.na(y))) stop("NA not permitted in response")

    nclass <- length(levels(y))
	
    ## Check for empty classes or small classes::
    if (any(table(y) == 0)) stop("Can't have empty classes in y.")
    #if (any(table(y) == 2)) stop("Can't have less than 3 cases per class in y.")

      moyenne <-double(nforest+1)   
      proba <-double(p)             
      maxiter=nforest
      nsample <- n               

      nrnodes <- 2 * trunc(nrow(x) / min(summary(y))) + 1


    ## Compiled code expects variables in rows and observations in columns.
    x <- t(x)
    storage.mode(x) <- "double"

    nt <- if (keep.forest) ntree else 1

        error.test <- double(1)
        rfout <- .C("classAgregTree",
                    x = x,
                    xdim = as.integer(c(p, n)),
                    y = as.integer(y),
                    nclass = as.integer(nclass),
                    Options = as.integer(c(
                    do.trace,
                    keep.inbag,
	            keep.forest, 
		    weight)),
                    ntree = as.integer(ntree),
		    nforest=as.integer(nforest),
		    maxiter=as.integer(maxiter), 
                    mtry = as.integer(mtry),
		    nstable=as.integer(nstable),      
                    nrnodes = as.integer(nrnodes),
                    ndbigtree = integer(ntree),
                    nodestatus = integer(nt * nrnodes),
                    bestvar = integer(nt * nrnodes),
                    treemap = integer(nt * 2 * nrnodes),
                    nodepred = integer(nt * nrnodes),
                    xbestsplit = double(nt * nrnodes),
                    inbag = if (keep.inbag) 
                    matrix(integer(n * ntree), n) else integer(n),
		    moyenne=moyenne,   
		    proba=proba,   
                    DUP=FALSE,
                    PACKAGE="ofw")[-1] 
        
	    if (keep.forest) {
            ## deal with the random forest outputs
            max.nodes <- max(rfout$ndbigtree)
            treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                             c(2, 1, 3))[1:max.nodes, , , drop=FALSE]
            }

	names(rfout$proba)=colnames(t(x))
	
	#compute the weights for each sample for output
	if(weight==T){
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

	
	
	
        cl <- match.call()
        cl[[1]] <- as.name("ofwCART")
        out <- list(call = cl,
                    type = "classification", 
                    classes = levels(y),
	            mean.error=matrix(rfout$moyenne[-1], ncol=1,dimnames=list(1:nforest, "MeanError")),
		    prob=rfout$proba,      #matrix(rfout$proba, ncol=1,dimnames=list(1:p, "probability")),
		    ntree = ntree,
		    nforest=nforest,
		    maxiter=rfout$maxiter, 
                    mtry = mtry,
	            do.trace=do.trace,
		    nstable=nstable,
		    weight=weight,
		    classWeight= if(!weight) NULL else classWeight,
		    sampleWeight= if(!weight) NULL else sampleWeight,
		    ##weightingOption= if (!weight) NULL else {list(classWeight=classWeight, sampleWeight=sampleWeight) },
                    forest = if (!keep.forest) NULL else {
                        list(ndbigtree = rfout$ndbigtree, 
                             nodestatus = matrix(rfout$nodestatus,
                             nc = ntree)[1:max.nodes,, drop=FALSE],
                             bestvar = matrix(rfout$bestvar, nc = ntree)[1:max.nodes,, drop=FALSE],
                             treemap = treemap,
                             nodepred = matrix(rfout$nodepred,
                             nc = ntree)[1:max.nodes,, drop=FALSE],
                             xbestsplit = matrix(rfout$xbestsplit,
                             nc = ntree)[1:max.nodes,, drop=FALSE],
                             nrnodes = max.nodes, ntree = ntree, nclass = nclass)
                    },
                    inbag = if (keep.inbag) rfout$inbag else NULL)
                   

    class(out) <- "ofwCART"
    return(out)
}
