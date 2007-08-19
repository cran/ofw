\name{ofwTune}
\alias{ofwTune}
\alias{ofwTune.default}


\title{Tuning the parameters for ofwCART or ofwSVM}
\description{
  \code{ofwTune} helps to tune the main parameters: \code{mtry}, \code{nforest} and \code{ntree} for ofwCART or \code{mtry} and \code{nsvm} for ofwSVM.
}

\usage{
\method{ofwTune}{default}(x, y, type="CART", ntree= if(type=="CART") 50 else NULL, 
	nforest= if(type=="CART") 100 else NULL, nsvm= if(type=="SVM") 500 else NULL, 
	mtry.test=seq(5,15,length=3), do.trace=FALSE, nstable=10, weight=FALSE, ...)
}
        
\arguments{
  \item{x}{A data frame with continuous values.}
  \item{y}{A response vector given as a factor (classification only).}
  \item{type}{The version of OFW to perform.}
  \item{ntree}{Number of trees to grow for each iteration (trees aggregation) in case ofwCART is chosen.}
  \item{nforest}{Total number of iterations to run if ofwCART.}
  \item{nsvm}{Total number of iterations to run if ofwSVM.}
  \item{mtry.test}{Vector defining the number of variables sampled according to the weight vector P to be tested. By default the vector is defined as 'seq(5,15,length=3)'.}
  \item{do.trace}{If set to some integer, then current iteration is printed for every \code{do.trace} iterations and the number of the first stable variables is output.}
  \item{nstable}{Need \code{do.trace} set to some integer. Stopping criterion before \code{nforest} or \code{nsvm} iterations are reached: if the \code{nstable} first weighted variables are the same after \code{do.stable} iterations, then stop. }
  \item{weight}{Should the weighting procedure be applied ?}
  \item{...}{not used currently.}
}

\details{
\code{ofwTune} consists in testing either ofwCART or ofwSVM with several variable subsets sizes in the given sequence \code{mtry.test}. For each \code{mtry}, the algorithm is performed twice and the function \code{ofwTune} outputs the intersection length of the first \code{nstable} variables selected with these 2 \code{ofw}. The value \code{mtry} to choose should be the one that gives the largest intersection. 

The total number of iteration to tune should then be 2 to 3 times the maximum iterations reached for the optimal \code{mtry}. 

In case of ofwCART, to choose \code{ntree}, the user should run \code{ofwTune} with several values of \code{ntree}. Usually, the more the trees the stabler the results.
}
 
\value{
A list with the following components:
  
  \item{type}{The classifier applied to ofw (either CART or SVM).}
  \item{nstable}{The number of stable variables chosen.}
  \item{mtry.test}{The different subset size tested.}
  \item{param}{A \code{2} by \code{length(mtry.test)} matrix indicating the number of stable variables obtained for each value of \code{mtry}.}
  \item{iter.max}{A \code{2} by \code{length(mtry.test)} matrix indicating the maximum number of iterations reached for each value of \code{mtry}.}
  \item{weight}{If \code{TRUE} the weighted procedure was performed during the learning.}
}

\note{
	The computation of \code{ofwTune} might be slow as it consists in launching ofw $2*$length(\code{mtry.Test}).% with an early stopping criterion (\code{nstable})

	%Usually, the optimal value of \code{mtry} should be the smallest that gives stable results.
	%\code{mtryTest} should not be greater than \code{nstable}.
	
	%\code{ntree=150} seems to give a good compromise between stable results and computation time.
	
}

\references{
	
  	L\^e Cao, K-A., Gon\c calves, O., Besse, P. and Gadat, S. (2007), \emph{Selection of biologically relevant genes with a wrapper stochastic algorithm} Statistical Applications in Genetics and Molecular Biology: Vol. 6: Iss.1, Article 29.
}

\author{Kim-Anh L\^e Cao \email{Kim-Anh.Le-Cao@toulouse.inra.fr}

	Patrick Chabrier \email{Patrick.Chabrier@toulouse.inra.fr}}

\seealso{\code{\link{ofw}}}

\examples{
## On data set "srbct"
#data(srbct)
#attach(srbct)
#tune.cart <- ofwTune(srbct, as.factor(class), type="CART", ntree=50, nforest=200, mtry.test=seq(5,10,length=2))
#tune.cart
#tune.svm <- ofwTune(srbct, as.factor(class), type="SVM", nsvm=500, mtry.test=seq(5,10,length=2))
#tune.svm
##Using do.trace options 
#tune.cart <- ofwTune(srbct, as.factor(class), type="CART", ntree=50, nforest=200, mtry.test=seq(5,10,length=2), do.trace=50, nstable=5)
#tune.cart
#tune.svm <- ofwTune(srbct, as.factor(class), type="SVM", nsvm=500, mtry.test=seq(5,10,length=2), do.trace=100, nstable=5)
#tune.svm

#detach(srbct)
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{tree}


