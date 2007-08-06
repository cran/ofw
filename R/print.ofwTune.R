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

"print.ofwTune" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("  Classifier applied to ofw: ", x$type, "\n",sep="")
  if(is.numeric(x$do.trace)) cat("  Stopping criterion was tested every ", x$do.trace," iterations \n",sep="") 

  if((is.numeric(x$do.trace)) &&  (x$itermax[1,which(x$param[2,]==max(x$param[2,]))[1]] < x$iteration)) {cat("  Each algorithm stopped when the first ", x$nstable," weighted variables were the same \n",sep="")} else {cat("  Each algorithm stopped when the total number of iteration ",x$iteration ," was reached\n",sep="")} 
  cat("  The optimal mtry can be ", x$param[1,which(x$param[2,]==max(x$param[2,]))[1]], "\n",sep="")
   cat("  Maximum iterations reached for this mtry value were ", x$itermax[1,which(x$param[2,]==max(x$param[2,]))[1]], " and ",x$itermax[2,which(x$param[2,]==max(x$param[2,]))[1]], " for both tries \n",sep="") 
  if(x$weight==FALSE ) {cat("  No weighted procedure \n")} else {cat("  Weighted procedure \n")}

}



