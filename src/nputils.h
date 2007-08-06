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

#ifndef NPU_H
#define NPU_H

#include "nptypes.h"

/* test if the bit at position pos is turned on */
#define isBitOn(x,pos) (((x) & (1 << (pos))) > 0)
/* swap two integers */
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))


int pack(int l, int *icat);
void unpack(int npack, int *icat);

void zeroInt(int *x, int length);
void zeroDouble(double *x, int length);
void prepare(int *cl, const int nsample, const int nclass, const int ipi, 
	     double *pi, double *pid, int *nc, double *wtt);



/** \brief constructs the matrix a an the matrix b
 * 
 * constructs the mdim by nsample integer array a.  For each 
 * numerical variable with values x(m, n), n=1, ...,nsample, the x-values 
 * are sorted from lowest to highest.  Denote these by xs(m, n).  Then 
 * a(m,n) is the case number in which xs(m, n) occurs. The b matrix is 
 * also contructed here.  If the mth variable is categorical, then 
 * a(m, n) is the category of the nth case number.
 * @param p_m_x blablabla.
 * @param p_m_a blablabla.
 * @param p_m_b blablabla.
 */
void makeA(p_matrixdouble_t p_m_x, 
	   p_matrixint_t p_m_a, 
	   p_matrixint_t p_m_b);

/** \brief blablabla
 * 
 * blablabla
 * @param p_m_a blablabla.
 * @param nuse blablabla.
 * @param ncase blabla.
 * @param jin blablabla.
 */
void modA(p_matrixint_t p_m_a , 
	  int *nuse, 
	  int *ncase, 
	  int *jin);

void Xtranslate(double *x, 
		int mdim, 
		int nrnodes, 
		int nsample, 
		int *bestvar, 
		int *bestsplit, 
		int *bestsplitnext,
		double *xbestsplit, 
		int *nodestatus, 
		int treeSize);

/* void permuteOOB(int m, double *x, int *in, int nsample, int mdim); */

/* void computeProximity(double *prox, int oobprox, int *node, int *inbag,  */
/*                       int *oobpair, int n); */


/* Node status */
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

#endif /* RF_H */
