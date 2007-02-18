\name{evaluate.learn}
\alias{evaluate}
\alias{evaluate.learn}


\title{Error rate assessment of ofw}
\description{
	The error rate assessment e.632+ (Efron and Tibshirani, 1997) is performed on ofw applied to \code{CART} or \code{SVM}. It requires first to launch \code{learn}. 
}

\usage{
\method{evaluate}{learn}(obj, maxvar=15, type=obj$type, nvar=if(obj$type=="CART")  obj$nclass+1 else NULL, ntreeTest= if(obj$type=="CART") 100 else NULL, weight=FALSE,...)
}
	    
	   
\arguments{
  \item{obj}{An object from class \code{learn}.}
  \item{maxvar}{Size of the evaluated variable selection.}
  \item{type}{Classifier used in the object from class \code{learn}}
  \item{nvar}{If \code{CART}, number of randomly sampled variables in the selection that are used to construct each tree. Should be at least \code{obj$nclass+1} to ensure generalizable trees.}
  \item{ntreeTest}{If \code{CART}, number of trees aggregated to evaluate the performance of ofwCART when each variable enters the selection.}
  \item{weight}{Should the weighting procedure be applied during the evaluation phase ?}

  \item{...}{not used currently.}
}

\details{
	In the case of data sets with a small number of samples (e.g microarray data), the use of e.632+ bootstrap error seems appropriate to assess the performance of the algorithm. With \code{CART}, as classification trees are unstable by nature, \code{ntreeTest} trees are aggregated. Furthermore, to avoid the wrong evaluation where the \code{ntreeTest} trees would be constructed only with the 'best' variables among the selected variables, \code{nvar} variables are randomly sampled in the evaluated selection to ensure a fair evaluation.
	}
	
\value{
  An object of class \code{evaluate}, which is a list with the
  following components:

  \item{maxvar}{Size of the evaluated variable selection.}
  \item{nvar}{Number of randomly sampled variables in the selection that are used to construct each tree.}
  \item{weight.eval}{Was the weighting procedure applied during the evaluation step ?}
  \item{weight.learn}{Was the weighting procedure applied during the learning step ?}
  \item{ntreeTest}{If \code{CART}, number of aggregated trees as variable enters the selection.}
  \item{matTrain}{A \code{nsample} by \code{Bsample} matrix indicating the training samples in each bootstrap sample.}
  \item{matProb}{A \code{nvariable} by \code{Bsample} matrix for each probability distribution learnt.}
  \item{error}{The evaluated e.632+ boostrap error as each variable enters the selection.}
  \item{sampleWeight}{if \code{weight==T}, the \code{n} by \code{Bsample} matrix indicating each sample weight in each bootstrap sample.} 
  \item{matPredInbag}{A \code{nvariable} by \code{Bsample} matrix indicating the prediction of the inbag samples.}
  \item{matPredTest}{A \code{nvariable} by \code{Bsample} matrix indicating the prediction of the test samples.}

}

\note{
	The e.632+ code comes from the \code{ipred} package.
	
	This type of evaluation should only be used to compare several methods and not to assess the performance of only one method.
	
	}

\references{
	Efron, B. and Tibshirani R.J. (1997), \emph{Improvements on cross-validation: the e.632+ bootstrap method}, Journal of American
Statistical Association 92, 548-560.

  	L\^e Cao, K-A., Gon\c calves, O., Besse, P. and Gadat, S. (2007), \emph{Selection of biologically relevant genes with a wrapper stochastic algorithm} Statistical Applications in Genetics and Molecular Biology: Vol. 6: Iss.1, Article 29.
}

\author{Kim-Anh L\^e Cao \email{Kim-Anh.Le-Cao@toulouse.inra.fr} \newline
Patrick Chabrier \email{Patrick.Chabrier@toulouse.inra.fr}}

\seealso{\code{\link{learn}}}

\examples{
## On data set "srbct"
#data(srbct)
#attach(srbct)
#learn.boot.cart <- learn(srbct, as.factor(class), type="CART", ntree=50, nforest=100, mtry=5, Bsample=3)
#eval.boot.cart <- evaluate(learn.boot.cart, ntreeTest=50, maxvar=10)
#plot(eval.boot.cart$error, type="l")
#learn.boot.svm <- learn(srbct, as.factor(class),type="SVM", nsvm=500, mtry=5, Bsample=3)
#eval.boot.svm <- evaluate(learn.boot.svm, maxvar=10)
#plot(eval.boot.svm$error, type="l")

#detach(srbct)
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{tree}


