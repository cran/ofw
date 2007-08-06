/*
 * Copyright (C) 2007 Kim-Anh Lê Cao, Patrick Chabrier, INRA,
 * French National Institute for Agricultural Research.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef CT_H
#define CT_H

#include "nptypes.h"


/** \brief constructs one classification tree
 *
 * @param a observation number when each variable is ordered
 * @param b factorization of the variable
 * @param class class of the observations
 * @param mDim number of variables
 * @param nsample number of observations
 * @param nClass number of classes
 * @param treemap is a nrnodes * 2 matrix indicating where the nodes were split cf predictClassTree
 * @param best=splitVar best variable on which to split
 * @param bestSplit observation where is the split that goes on the left
 * @param bestSplitNext next observation where is the split
 * @param tgini decrease Gini on the variable that is split
 * @param nodeStatus status of the ndoe (1, 2, or -1)
 * @param nodePop observations number that are in each node
 * @param nodeStart node number on which to start
 * @param tclassCount matrix nClass*nrnodes count number of observations of each class for each node
  * @param classCount same as tclassCount but as a vector (one node)
 * @param ta cf movedata is a matrix transfer (besoin ici?)
 * @param nrnodes maximum number of nodes (fixed)
 * @param idmove cf movedata Pas besoin ici
 * @param ndsize minimum number of observation in a node (fixed value)
 * @param ncase given by movedata
 * @param mtry number of variables randomly drawn to choose the split
 * @param varUsed if 1 this variable has been used for splitting a node
 * @param nodeclass assigned class of the node
 * @param ndbigtree number of nodes in the tree
 * @param weight observation weight
 * @param wr number of observation per class on the right node (weighted)
 * @param wl  number of observation per class on the left node (weighted)
 * @param nuse computed in modA number of observations in bag
 */

void buildtree(int *a, 
	       int *b, 
	       int *class, 
	       int *mDim, 
	       int *nsample,
	       int *nClass, 
	       int *treemap,  
	       int *best,  
	       int *bestSplit, 
	       int *bestSplitNext,
	       double *tgini, 
	       int *nodeStatus, 
	       int *nodePop, 
	       int *nodeStart,  
	       double *tclassCount,//tableau
	       double *classCount,//vecteur
	       int *ta, //utilise ds movedata
	       int *nrnodes, //nb noeud max ds arbre,fixé par defaut
	       int *idmove, //indic d'obs a gche ou drte
	       int *ndsize, //val donnee par defaut
	       int *ncase,
	       int *mtry,
	       int *varUsed,  //equiv de iv
	       int *nodeclass, //classe assignee au noeud
	       int *ndbigtree, //nb noeuds existants ds larbre
	       double *weight, //equiv de win
	       double *wr, 
	       double *wl,
	       int *nuse, //valeur par defaut
	       int *omega
	       );

/** \brief reorder a once the split is chosen and performed
 *
 * @param a observation number when each variable is ordered
 * @param ta is a transfer matrix to modify a
 * @param mDim number of variables
 * @param nsample number of observations
 * @param ndStart index where to begin (cf builtree.c)
 * @param ndEnd index where to finish
 * @param idmove 1 if observation goes on the left, 0 otherwise
 * @param ncase reorder the observations after the split is done
 * @param splitVar variable on which to split
 * @param nbest observation where to split
 * @param ndendl = nbest utile ici?
 * @return ncase ndendl and a
 */

void movedata(int *a, 
	      int *ta,
	      int *mDim, 
	      int *nsample, 
	      int *ndStart, 
	      int *ndEnd,
	      int *idmove, 
	      int *ncase, 
	      int *splitVar, 
	      int *nbest,
	      int *ndendl);

/** \brief finds the best variable on which to split the node in the tree
 *
 * @param a observation number when each variable is ordered
 * @param b factorization of the variable
 * @param class class of the observations
 * @param mDim number of variables
 * @param nClass number of classes
 * @param ndStart index where to begin (cf builtree.c)
 * @param ndEnd index where to finish
 * @param classCount number of observations for each class in a node
 * @param splitVar variable on which to split
 * @param decGini decrease Gini
 * @param nbest observation where to split
 * @param splitStatus node status
 * @param mtry number of variables randomly drawn to choose the split
 * @param weight observation weight
 * @param wr number of observation per class on the left node (weighted)
 * @param wl number of observation per class on the right node (weighted)
 * @return splitVar decGini nbest and splitStatus (verifier)
 */
 
void findbestsplit(int *a, 
		   int *b, 
		   int *class, 
		   int *mDim, 
		   int *nClass, 
		   int *ndStart, 
		   int *ndEnd, 
		   double *classCount,
		   int *splitVar, 
		   double *decGini, 
		   int *nbest, 
		   int *splitStatus, 
		   int *mtry, /* racine carrée du nombre de variables */
		   double *weight, 
		   double *wr, 
		   double *wl,
		   int *omega
		   );




#endif
