\name{learn}
\alias{learn}
\alias{learn.default}


\title{Learning ofw for error rate assessement}
\description{
  \code{learn} simply performs the learning step of ofw on several bootstrap samples in order to assess the error rate on the test observations (see \code{evaluate})}
}

\usage{
\method{learn}{default}(x, y, type="CART", ntree= if(type=="CART") 50 else NULL, nforest= if(type=="CART") 100 else NULL, niteration= if(type=="SVM") 20000 else NULL, mtry=5, do.trace=FALSE, nstable=50, weight=FALSE, Bsample=5, ...)
}
        
\arguments{
  \item{x}{A data frame with continuous values.}
  \item{y}{A response vector given as a factor (classification only).}
  \item{type}{Classifier used: either \code{CART} or \code{SVM}.}
  \item{ntree}{If \code{CART}, number of trees to grow for each iteration (trees aggregation).}
  \item{nforest}{If \code{CART}, number of iterations to run. This should not be set to too small a number, to ensure the convergence of the algorithm.}
  \item{niteration}{If \code{SVM}, number of iterations to run. This should be set to a very large number, to ensure the convergence of the algorithm.}
  \item{mtry}{Number of variables sampled according to the weight vector P as candidates for each tree or SVM. This should be small enough to ensure stable results of the algorithm.}
  \item{do.trace}{If set to some integer, then current iteration is printed for every \code{do.trace} iterations and the number of the first stable variables is output.}
  \item{nstable}{Need \code{do.trace} set to some integer. Stopping criterion before \code{nforest} iterations are reached: if the \code{nstable} first weighted variables are the same after \code{do.stable} iterations, then stop. }
  \item{weight}{Should the weighting procedure be applied ?}
  \item{Bsample}{Number of bootstrap samples for the learning step.}
  
  \item{...}{not used currently.}
}

\details{The object from class \code{learn} will be used in the generic function \code{evaluate}}
 
\value{
  An object of class \code{learn}, which is a list with the
  following components:

  \item{x}{Original input}
  \item{y}{Original predictor}
  \item{type}{Classifier used: either \code{CART} or \code{SVM}}
  \item{nsample}{Total number of samples}
  \item{nclass}{Number of levels of \code{y} (number of classes)}
  \item{weight.learn}{If \code{TRUE} the weighted procedure was performed during the learning.}
  \item{Bsample}{Number of bootstrap samples on which ofwCART is learnt.}
  \item{matTrain}{A \code{n} by \code{Bsample} matrix indicating the training samples in each bootstrap sample.}
  \item{matProb}{A \code{nvariable} by \code{Bsample} matrix for each probability distribution learnt.}
  \item{weightingOption}{if \code{weight==T} : \newline
  	-classWeight: a \code{nclass} by \code{Bsample} matrix indicating classes weights of each bootstrap sample.\newline
  	-sampleWeight: a \code{nsample} by \code{Bsample} matrix indicating each sample weight in each bootstrap sample.
 	}
}

\note{

	The computation of learn is slow as it requires to launch \code{ofw} \code{Bsample} times.

	Parallelized computations are possible with \code{ofw} and  \code{Rmpi} library for the learning step and \code{evaluateCARTparallel} for the evaluation step.
}

\references{
	
  	L\^e Cao, K-A., Gon\c calves, O., Besse, P. and Gadat, S. (2007), \emph{Selection of biologically relevant genes with a wrapper stochastic algorithm} Statistical Applications in Genetics and Molecular Biology: Vol. 6: Iss.1, Article 29.
  	L\^e Cao, K-A., Bonnet, A., and Gadat, S., \emph{Multiclass classification  and gene selection with a stochastic algorithm}  \url{htp://www.lsp.ups-tlse.fr/}, Publication LSP\_2007\_02.
}

\author{Kim-Anh L\^e Cao \email{Kim-Anh.Le-Cao@toulouse.inra.fr}\newline
Patrick Chabrier \email{Patrick.Chabrier@toulouse.inra.fr}}

\seealso{\code{\link{evaluate.learn}}}

\examples{
## On data set "data"
#data(srbct)
#attach(srbct)
#learn.boot.cart <- learnCART(srbct, as.factor(class),type="CART", ntree=50, nforest=200, mtry=5, Bsample=3)
#learn.boot.svm <- learnSVM(srbct, as.factor(class), type="SVM", niteration=500, mtry=5, Bsample=3)
#detach(srbct)
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{tree}


