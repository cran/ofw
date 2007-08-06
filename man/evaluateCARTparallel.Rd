\name{evaluateCARTparallel}
\alias{evaluateCARTparallel}



\title{Error rate assessment of ofwCART}
\description{
	The error rate assessment e.632+ (Efron and Tibshirani, 1997) is performed on ofwCART. This second version of \code{evaluateCART} allows to perform the learning step independantly (for example with parallel computing).
}

\usage{
evaluateCARTparallel(x, y, matTrain, matProb, maxvar=15, nvar=nlevels(y)+1,
	 ntreeTest=100, weight=F)
}

	   
\arguments{
  \item{x}{A data frame with continuous values.}
  \item{y}{A response vector given as a factor (classification only).}
  \item{matTrain}{A \code{nsample} by \code{Bsample} matrix indicating the training samples in each bootstrap sample.}
  \item{matProb}{A \code{nvariable} by \code{Bsample} matrix for each probability distribution learnt.}
  \item{maxvar}{Size of the evaluated variable selection.}
  \item{nvar}{Number of randomly sampled variables in the selection that are used to construct each tree. Should be at least \code{nlevels(y)+1} to ensure generalizable trees.}
  \item{ntreeTest}{Number of trees aggregated to evaluate the performance of ofwCART when each variable enters the selection.}
  \item{weight}{Should the weighting procedure be applied during the evaluation phase ?}

%  \item{...}{optional parameters to be passed to the low level function \code{}.}  %enlever
}

\details{
	see \code{evaluateCART}
	}
	

\value{
  An object of class \code{evaluateCARTparallel}, which is a list with the
  following components:

  \item{maxvar}{Size of the evaluated variable selection.}
  \item{nvar}{Number of randomly sampled variables in the selection that are used to construct each tree.}
  \item{weight.eval}{Was the weighting procedure applied during the evaluation step ?}
  \item{ntreeTest}{Number of aggregated trees as each variable enters the selection.}
  \item{matTrain}{A \code{nsample} by \code{Bsample} matrix indicating the training samples in each bootstrap sample.}
  \item{matProb}{A \code{nvariable} by \code{Bsample} matrix for each probability distribution learnt.}
  \item{error}{The evaluated e.632+ boostrap error as each variable enters the selection.}
  \item{sampleWeight}{if \code{weight==T}, a \code{n} by \code{Bsample} matrix indicating each sample weight in each bootstrap sample.} 
  \item{matPredInbag}{A \code{nvariable} by \code{Bsample} matrix indicating the prediction of the inbag samples.}
  \item{matPredTest}{A \code{nvariable} by \code{Bsample} matrix indicating the prediction of the test samples.}

   

}

\note{
	The e.632+ code comes from the \code{ipred} package from Thorsten.
	
	This type of evaluation should only be used to compare several methods and not to assess the performance of only one method.
	
	}

\references{
	Efron, B. and Tibshirani R.J. (1997), \emph{Improvements on cross-validation: the e.632+ bootstrap method}, Journal of American
Statistical Association 92, 548-560.

  	L\^e Cao, K-A., Gon\c calves, O., Besse, P. and Gadat, S. (2007), \emph{Selection of biologically relevant genes with a wrapper stochastic algorithm} Statistical Applications in Genetics and Molecular Biology: Vol. 6: Iss.1, Article 29.
}

\author{Kim-Anh L\^e Cao \email{Kim-Anh.Le-Cao@toulouse.inra.fr} \newline
Patrick Chabrier \email{Patrick.Chabrier@toulouse.inra.fr}.}

\seealso{\code{evaluate}}

\examples{
## On data set "data"
##first learn ofwCART on Bsample bootstrap samples and store matTrain.data and matProb.data
##data.evalCARTparallel <- evaluateCARTparalell(data[,-1], as.factor(data[,1], matTrain=matTrain.data, matProb=matProb.data, ntreeTest=50, maxvar=10)
##plot(data.evalCARTparallel$error, type="l")
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{tree}


