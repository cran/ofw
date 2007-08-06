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

getTree <- function(ofwobj, k=1) {
  if (is.null(ofwobj$forest)) {
    stop("No forest component in ", deparse(substitute(ofwobj)))
  }
  if (k > ofwobj$ntree) {
    stop("There are fewer than ", k, "trees in the last forest")
  }
      tree <- cbind(ofwobj$forest$treemap[,,k],
                    ofwobj$forest$rightDaughter[,k],
                    ofwobj$forest$bestvar[,k],
                    ofwobj$forest$xbestsplit[,k],
                    ofwobj$forest$nodestatus[,k],
                    ofwobj$forest$nodepred[,k])[1:ofwobj$forest$ndbigtree[k],]


  dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
                                         "split var", "split point",
                                         "status", "prediction"))

  tree
}

