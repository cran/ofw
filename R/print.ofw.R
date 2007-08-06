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

"print.ofw" <-
function(x, ...) {

if(x$type=="CART"){
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("  Classifier: ", x$type, "\n",sep="")
  #cat("  Number of trees: ", x$ntree, "\n",sep="")
  if(is.numeric(x$do.trace) &&  (x$maxiter < x$niteration)) {cat("  Number of iterations reached: ", x$maxiter, "\n",sep="")} else {cat("  Number of iterations: ", x$nforest, "\n",sep="")}
  
  cat("  No. of variables tried at each iteration: ", x$mtry, "\n\n", sep="")
  if(x$weight==FALSE ) {cat("  No weighted procedure \n")} else {cat("  Weighted procedure \n")}


  if(is.numeric(x$do.trace)) cat("  Stopping criterion was tested every ", x$do.trace," iterations \n",sep="") 
  if((is.numeric(x$do.trace)) && (x$maxiter < x$niteration)) cat("  The algorithm stopped when the first ", x$nstable," weighted variables were the same \n",sep="")
} 


if(x$type=="SVM"){
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("  Classifier: ", x$type, "\n",sep="")
  #cat("  Number of trees: ", x$ntree, "\n",sep="")
  if(is.numeric(x$do.trace) &&  (x$maxiter < x$niteration)) {cat("  Number of iterations reached: ", x$maxiter, "\n",sep="")} else {cat("  Number of iterations: ", x$niteration, "\n",sep="")}
  
  cat("  No. of variables tried at each iteration: ", x$mtry, "\n\n", sep="")
  if(x$weight==FALSE ) {cat("  No weighted procedure \n")} else {cat("  Weighted procedure \n")}


  if(is.numeric(x$do.trace)) cat("  Stopping criterion was tested every ", x$do.trace," iterations \n",sep="") 
  if((is.numeric(x$do.trace)) && (x$maxiter < x$niteration)) cat("  The algorithm stopped when the first ", x$nstable," weighted variables were the same \n",sep="")
}

}
