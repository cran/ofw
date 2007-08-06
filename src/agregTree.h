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

#ifndef AGREGTREE_H
#define AGREGTREE_H



/** \brief Computes the out of bag error rate
 * 
 * blablabla
 * @param nsample number of observations
 * @param nclass number of classes
 * @param cl class of the observations
 * @param jerr if 1 indicates wrong predicted class for each observation
 * @param counttr counts the number of votes from each tree for each observation
 * @param out indicates if observation is in the out bag
 * @param errtr number of prediction error for each observation
 * @param jest predicted class (?)
 * @param cutoff cutoff value for the predicted class
 * @return jerr jest and errtr
 */

//void oob(int nsample, int nclass,  int *cl,int *jerr,
//         int *counttr, int *out, double *errtr, int *jest, double *cutoff);




/** \brief Computes the out of bag error rate errtr for each tree
 * 
 * @param nsample number of observations
 * @param cl class of the observations
 * @param jtr the predicted class of the oob cases
 * @param jin indicates if observation is in the out bag
 * @param errtr number of prediction error for each observation
 * @return errtr
 */

void oobTree(
	 int nsample, 
	 int *cl, 
	 int *jtr, 
	 int *jin, 
	 double *errtr, //taux d erreur (sortie),
	 int jb
	 ) ;


/** \brief Computes the weighted out of bag error rate errtr for each tree
 * 
 * @param nsample number of observations
 * @param cl class of the observations
 * @param jtr the predicted class of the oob cases
 * @param jin indicates if observation is in the out bag
 * @param errtr number of prediction error for each observation
 * @param jb the tree number in the forest loop
 * @param classWeight the weight on each class 
 * @return errtr
 */
void oobTreeWeight(
	 int nsample, 
	 int *cl, 
	 int *jtr,  
	 int *jin, 
	 double *errtr, 
	 int jb,
	 double *classWeight
	 );


/** \brief predict the class of each observation in the tree
 * 
 * @param x input data
 * @param n number of observations
 * @param mdim number of variables
 * @param treemap is a nrnodes * 2 matrix indicating where the nodes were split
 * @param nodestatus status of the node 
 * @param xbestsplit cf nputils.c the best observation where to split
 * @param bestvar the best variable on which to split
 * @param nodeclass assigned class of the terminal node
 * @param jts predicted class of the observations
 * @param nodex node number where each observation lands
 * @return nodex and jts
 */
void predictClassTree(double *x, 
		      int n, 
		      int mdim, 
		      int *treemap,
		      int *nodestatus, 
		      double *xbestsplit,
		      int *bestvar, 
		      int *nodeclass,
		      int *jts, 
		      int *nodex);


/** \brief Construction of the forest
 * 
 */
void classAgregTree(double *x, 
	     int *dimx, 
	     int *cl, 
	     int *ncl, 
	     int *Options, 
	     int *ntree, 
	     int *nforest, 
	     int *maxiter,
	     int *nvar,
	     int *nstable,
	     int *nrnodes, 
	     int *ndbigtree, 
	     int *nodestatus, 
	     int *bestvar, 
	     int *treemap, 
	     int *nodeclass, 
	     double *xbestsplit,
             int *inbag, 
	     double *moyenne, 
	     double *proba
	     );

void classForest(p_matrixdouble_t p_m_x, 
		 int *dimx, 
		 int *cl, 
		 int *ncl, 
		 int *Options, 
		 int *ntree, 
		 int *nforest,
		 int *nvar,
		 double *classwt, 
		 int *nrnodes, 
		 int *ndbigtree, 
		 int *nodestatus, 
		 int *bestvar, 
		 int *treemap, 
		 int *nodeclass, 
		 double *xbestsplit, 
		 double *errtr, 
		 int *inbag,
		 double *tgini,
		 double *wl,
		 double *wr,
		 double *classpop,   
		 double *tclasspop, 
		 double *tx, 
		 double *win,
		 double *tp,
		 int * out,
		 int * bestsplitnext,
		 int * bestsplit,
		 int * nodepop,
		 int * nodestart,
		 int * jin,
		 int * nodex,
		 int * ta,
		 int * ncase,         
		 int * jerr,
		 int * varUsed,
		 int * jtr,
		 int * classFreq ,
		 int * idmove,     
		 p_matrixint_t p_m_at,
		 p_matrixint_t p_m_a,
		 p_matrixint_t p_m_b,
		 double * proba,
		 int *omega,
		 double *gradient,
		 int iter,
		 double * moyenne,
		 double * sumProba,
		 int *perm,
		 int * inside,
		 int * numWeight,
		 double * classWeight
		 );

/** \brief initialize proba as a uniform distribution on mDim
 * 
 * @param proba the probability on mDim initialized as a uniform distribution
 * @param mDim number of variables
 * @return uniform distribution on mDim
 */
void proba_init(double *proba, int *mDim);


/** \brief returns the variable subset omega drawn with resp. to proba and of size s_omega
 * 
 * @param proba probability on the variables
 * @param mdim number of variables
 * @param s_omega size of the subset of variables drawn with respect to proba
 * @return omega: the subset of variables drawn with resp. to proba
 */
void draw_omega(double *proba,
		int *mDim,
		int *omega,
		int *s_omega,
		double * sumProba,
		int *perm
 		);


/** \brief project the probability on the simplex 
 * 
 * @param proba probability on the variables
 * @param mdim number of variables
 * @return probability projected on the simplex
 */
void simplexe( double *proba,
	       int *mDim);


#endif
