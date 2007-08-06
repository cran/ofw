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

"ofwCART.formula" <-
    function(formula, data = NULL, ..., subset, na.action = na.fail) {	
### formula interface for randomForest.
### code gratefully stolen from svm.formula (package e1071).
### CHANGER POUR randomForest.formula ??



    if (!inherits(formula, "formula"))
        stop("method is only for formula objects")
    m <- match.call(expand = FALSE)
    names(m)[2] <- "formula"
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    y <- model.response(m)
    if(!is.null(y)) m <- m[, -1, drop=FALSE]
    for (i in seq(along=ncol(m))) {
        if(is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
    }
    ret <- ofwCART(m, y, ...)
    cl <- match.call()
    cl[[1]] <- as.name("ofwCART")
    ret$call <- cl
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action"))) 
        ret$na.action <- attr(m, "na.action")
    class(ret) <- c("ofwCART.formula", "ofwCART")
    return(ret)
}
