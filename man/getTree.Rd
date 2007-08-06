\name{getTree}
\alias{getTree}
\title{Extract a single tree from the last iteration of ofwCART.}
\description{
  This function extract the structure of a tree from a
  \code{ofwCART} object for the last iteration. 
}
\usage{
getTree(ofwobj, k=1)
}
\arguments{
  \item{ofwobj}{a \code{\link{ofw}} object applied with classifier CART.}
  \item{k}{which tree to extract form the last forest?}
}
\value{
  A matrix (or data frame, if \code{labelVar=TRUE}) with six columns and
  number of rows equal to total number of nodes in the tree.  The six
  columns are: 
  \item{left daughter}{the row where the left daughter node is; 0 if the
  node is terminal}
  \item{right daughter}{the row where the right daughter node is; 0 if
    the node is terminal}
  \item{split var}{which variable was used to split the node; 0 if the
    node is terminal}
  \item{split point}{where the best split is; see Details for
    categorical predictor}
  \item{status}{is the node terminal (-1) or not (1)}
  \item{prediction}{the prediction for the node; 0 if the node is not
    terminal}
}
\details{
 This function comes from the \code{randomForest} package.
 }
%\references{
%}
\seealso{
  \code{\link{ofw}}
}
\examples{
data(srbct)
attach(srbct)
##ofwCART
learn.cart.keep <- ofw(srbct, as.factor(class),type="CART", ntree=50, nforest=10, mtry=5, keep.forest=TRUE)
getTree(learn.cart.keep, k=3)
detach(srbct)
}
\author{}
\keyword{tree}
