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


#this function learns the SVM on a bootstrap sample and computes the OOB error rate


svmlearn=function(x, y,lw, P, nvariable, nsample, weight, vectWeight = NULL)
{
	G=vector(length=length(P))
	W=sample(nvariable, size=lw, prob=P, replace=FALSE) #tirer W selon la proba P

	cont=TRUE
	while (cont){
		train=sample(1: nsample, nsample, replace=TRUE)
		if ((any(table(y[train]) == 0))|| (any(table(y[setdiff(1:nsample, train)]) == 0))) {cont=TRUE} else {cont=FALSE}
		}

	test = setdiff(1:nsample, train)                 #echantillon oob     
	data.train = x[train,W]       
	data.test=x[test,W]
	svm.train=svm(data.train, y[train], kernel="linear")
	mat=table(y[test],predict(svm.train, data.test))
	if (weight) {erreur=sum((apply(mat, 1, sum) -diag(mat))*vectWeight)/sum(apply(mat, 1, sum)*vectWeight)} else {erreur= (length(test) - sum(diag(mat)))/length(test)}
	G[W]= G[W] +erreur/(1000*P[W])  
	return(list(G,erreur))
}



