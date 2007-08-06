\name{ofw}
\alias{ofw}
\alias{ofw.default}
\alias{print.ofw}
%\alias{ofwCART}
%\alias{ofwCART.default}
%\alias{print.ofwCART}

\title{Optimal Feature Weighting with CART and SVM}
\description{
  \code{ofw} implements the Optimal Feature Weighting stochastic algorithm for multiclass classification with Classification And Regression Trees (CART) or Support Vector Machine (SVM) as a classifier, in the context of continuous variables.
}
\usage{
\method{ofw}{default}(x,y,type="CART", ntree= if(type=="CART") 50 else NULL, nforest= if(type=="CART") 100 else NULL, niteration= if(type=="SVM") 100 else NULL, mtry=5, do.trace=FALSE, nstable=25, keep.inbag=if(type=="CART")  FALSE else NULL, keep.forest=if(type=="CART") TRUE else NULL, weight=FALSE, ...)
\method{print}{ofw}(x, ...)
}
        
\arguments{
  \item{x}{A data frame with continuous values (for the \code{print} method, an \code{ofw} object).}
  \item{y}{A response vector given as a factor (classification only).}
  \item{type}{Classifier used: either \code{CART} or \code{SVM}}
  \item{ntree}{If \code{CART}, number of trees to grow for each iteration (trees aggregation).}
  \item{nforest}{If \code{CART}, number of iterations to run. This should not be set to too small a number, to ensure the convergence of the algorithm.}
  \item{niteration}{If \code{SVM}, number of iterations to run. This should be set to a very large number, to ensure the convergence of the algorithm.}
  \item{mtry}{Number of variables sampled according to the weight vector P as candidates for each tree or SVM. This should be small enough to ensure stable results of the algorithm.}
  \item{do.trace}{If set to some integer, then current iteration is printed for every \code{do.trace} iterations and the number of the first stable variables is output.}
  \item{nstable}{Need \code{do.trace} set to some integer. Stopping criterion before \code{nforest} or \code{niterations} iterations are reached: if the \code{nstable} first weighted variables are the same after \code{do.stable} iterations, then stop. }
  \item{keep.inbag}{If \code{CART}, should an \code{n} by \code{ntree} matrix be returned that keeps track of which samples are ``in-bag'' in which trees (and how many times as it is a sampling with replacement) in the last forest.}
  \item{keep.forest}{If \code{CART}, and if set to \code{TRUE}, the last forest (or last iteration) will be retained in the output object and the \code{getTree} function can be used to see how the trees were constructed.}
  \item{weight}{Should the weighting procedure be applied ?}
  
  \item{...}{not used currently.}
}

\details{The Optimal Feature Weighting algorithm consists in learning the probability distribution P on all variables. The more useful the variables in the classification task, the heavier their weight. When the CART classifier is used, the trees are aggregated at each iteration (bagging).}
                   
\value{
  An object of class \code{ofw}, which is a list with the
  following components:

  \item{call}{The original call to \code{ofwCART}.}
  \item{type}{Classifier used.}
  \item{classes}{Level attributes of the classes in the \code{y} predictor.}
  \item{mean.error}{Internal mean error rate vector of length \code{nforest}.}
  \item{prob}{Probability distribution vector P or weighting vector of length the total number of variables. }
  \item{ntree}{If \code{CART}, number of trees grown for each iteration.}
  \item{nforest}{If \code{CART}, number of iterations asked in the procedure.}
  \item{niteration}{If \code{SVM}, number of iterations asked in the procedure.}
  \item{maxiter}{Actual number of iterations performed (different than \code{nforest} or \code{niteration} if \code{nstable} variables are obtained after \code{do.trace} steps).}
  \item{mtry}{Number of predictors sampled with P for each tree for each iteration.}
  \item{do.trace}{The number of iteration is printed every \code{do.trace} and the stopping criterion is tested.}
  \item{nstable}{Number of stable variables obtained if \code{maxiter < nforest}.}
  \item{weight}{If \code{TRUE} the weighted procedure was performed.}
  \item{weightingOption}{If \code{weight==T} : \newline
  	-classWeight: class weight vector.\newline
  	-sampleWeight: sample weight vector.}
  \item{forest}{If \code{CART}, a list that contains the entire last forest of the last iteration; \code{NULL} if \code{keep.forest=FALSE}.}
  \item{inbag}{If \code{keep.inbag=TRUE} a \code{n} by \code{ntree} matrix is returned with ``in-bag'' samples in which trees (and sampled how many times).}

  	
}

\note{
	The \code{ofw} for \code{CART} structure has been first largely inspired from the \code{randomForest} package which was itself based on Leo Breiman and Adele Cutler's Fortran code. The code now written in C (or R).
	
	The actual implementation of ofw is restrained to classification task with continuous variables. It has been especially developped to deal with p$>>$n data sets, such as microarrays.
	
	Normalisation has first to be performed by the user. For extremely large data sets, a large pre processing is advised to speed up the procedure. 
	
	In contrary to \code{CART}, the ofw version with \code{SVM} does not provide the internal mean error rate.

}

\references{
	Gadat, S. and Younes, L. (2007), \emph{A Stochastic Algorithm for Feature Selection in Pattern Recognition}, Journal of Machine Learning 8, 509-548
	
  	L\^e Cao, K-A., Gon\c calves, O., Besse, P. and Gadat, S. (2007), \emph{Selection of biologically relevant genes with a wrapper stochastic algorithm} Statistical Applications in Genetics and Molecular Biology: Vol. 6: Iss.1, Article 29.
  	
  	L\^e Cao, K-A., Bonnet, A., and Gadat, S., \emph{Multiclass classification  and gene selection with a stochastic algorithm}  \url{htp://www.lsp.ups-tlse.fr/}, Publication LSP\_2007\_02.
}

\author{Kim-Anh L\^e Cao \email{Kim-Anh.Le-Cao@toulouse.inra.fr} \newline
	Patrick Chabrier \email{Patrick.Chabrier@toulouse.inra.fr} \newline
	Based on Leo Breiman and Adele Cutler Fortran code and \code{randomForest} package for \code{ofwCART} version.}

\seealso{\code{\link{getTree}}, \code{\link{learn}}, \code{\link{evaluate.learn}}}

\examples{
## On data set "srbct"
data(srbct)
attach(srbct)
##ofwCART
learn.cart <- ofw(srbct, as.factor(class),type="CART", ntree=100, nforest=200, mtry=5)
print(learn.cart)
## Look at variable importance:
learn.cart$prob
## Look if the internal mean error is decreasing w.r.t the number of iterations:
plot(learn.cart$mean.error, type="l")
##ofwSVM
learn.svm <- ofwSVM(srbct, as.factor(class),type="SVM", niteration=500, mtry=5)
print(learn.svm)
## Look at variable importance:
learn.svm$prob
detach(srbct)
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{tree}
