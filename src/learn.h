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

#include "nptypes.h"
#include "nputils.h"


#ifndef E_H
#define E_H

/** \brief Construction of the forest
 * 
 */
void classLearn(double *x, 
		int *dimx, 
		int *cl, 
		int *ncl, 
		int *ntree,
		int *nvar,  //nb de variables tirées
		int *ngene,  //nb de var evaluees
		double *classwt, 
		double *cut, 
		int *counttr, 
		int *nrnodes, 
		int *ndbigtree, 
		int *nodestatus, 
		int *bestvar, 
		int *treemap, 
		int *nodeclass, 
		double *xbestsplit, 
		double *xts, 
		int *clts, 
		int *nts, 
	     double *countts,
		int *outclts, 
		int *outcl,
		double *errts,
		double *errin,
		double *proba,
		int *weight,
		double *classWeight
		);



/** \brief Computes the weighted or not  out of bag error rate errtr for each tree
 * 
 * @param nsample number of observations
 * @param nclass number of classes
 * @param cl class of the observations
 * @param p_m_counttr counts the number of times each tree votes for a predicted class for each observation
 * @param errtr number of prediction error for each observation
 * @param jest the predicted class of the oob cases
 * @param cutoff is the value for majority vote
 * @param jtr is the class node where each obs lands
 * @param nvote=jb is the number of trees that vote
 * @param weight if 1 the error rate is weighted
 * @param sampleWeight the weight on each sample 
 * @return errtr and jest
 */

void inbagError(
	 int nsample, 
	 int nclass,
	 int *cl, 
	 p_matrixint_t p_m_counttr,  //int *counttr,   
	 double *errtr, //taux d erreur (sortie)
	 int *jest,  //equiv outcl: classe predite
	 double *cutoff,
	 int *jtr,
	 int nvote,
	 int weight,
	 double *classWeight 
	 );

/** \brief Computes the test set error rate 
 * Probablement a modifier
 */
void TestSetError(p_matrixdouble_t p_m_countts, int *jts, int *clts, int *jet,
		  int ntest,int nclass, int nvote, double *errts,  
		   double *cutoff,
		  int weight,
		  double *classWeight );


#endif
