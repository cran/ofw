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

#include <R.h>
#include "nputils.h"

void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int)); 
}

void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double)); 
}



void makeA(p_matrixdouble_t p_m_x,
	   p_matrixint_t p_m_a, 
	   p_matrixint_t p_m_b) {
    /* makeA() constructs the mdim by nsample integer array a.  For each 
       numerical variable with values x(m, n), n=1, ...,nsample, the x-values 
       are sorted from lowest to highest.  Denote these by xs(m, n).  Then 
       a(m,n) is the case number in which xs(m, n) occurs. The b matrix is 
       also contructed here.  If the mth variable is categorical, then 
       a(m, n) is the category of the nth case number. */

    /* "traduction" a contient (par colonne) l'indice de la valeur
       de x dans le tableau ordonnée
    */

    /* alors que dans la matrice b on substitue les valeurs pas des indices
       de classes */
    int i, j, n1, n2, *index;

    int nsample = p_m_x->nrows;
    int mdim = p_m_x->ncols;
      
    double *v;

    v     = (double *) Calloc(nsample, double);
    index = (int *) Calloc(nsample, int);
    
    for (i = 0; i < mdim; ++i) {
	for (j = 0; j < nsample; ++j) {
	  v[j] = p_m_x->array[j][i];
	  index[j] = j + 1;
	}
	
	/*  this sorts the v(n) in ascending order. index(n) is the case 
	    number of that v(n) nth from the lowest (assume the original 
	    case numbers are 1,2,...).  */
	R_qsort_I(v, index, 1, nsample);
	
	for (j = 0; j < nsample-1; ++j) {
	  n1 = index[j];
	  n2 = index[j + 1];
	  p_m_a->array[j][i] = n1;
	  
	  /* first class is 1 */
	  if (j == 0)  p_m_b->array[n1-1][i] = 1;
	  
	  /* either the class change or not */
	  p_m_b->array[n2-1][i] =  (v[j] < v[j + 1]) ?
	    p_m_b->array[n1-1][i] + 1 : p_m_b->array[n1-1][i];
	}
	
	p_m_a->array[nsample-1][i] = index[nsample-1];
    }
    Free(index);
    Free(v);
}


void modA(p_matrixint_t p_m_a, 
	  int *nuse,
	  int *ncase, 
	  int *jin) {
  
  int i, j, k, m, nt;
  
  int nsample = p_m_a->nrows;
  int mdim = p_m_a->ncols;
  
  *nuse = 0;
  for (i = 0; i < nsample; ++i) if (jin[i]) (*nuse)++;
  
  for (i = 0; i < mdim; ++i) {
    k = 0;
    nt = 0;
    for (j = 0; j < nsample; ++j) {
      if (jin[p_m_a->array[k][i] - 1]) {
	p_m_a->array[nt][i] = p_m_a->array[k][i];
	k++;
	} else {
	for (m = 0; m < nsample - k; ++m) {
	    if (jin[p_m_a->array[k + m][i] - 1]) {
	      p_m_a->array[nt][i] = p_m_a->array[k + m][i];
	      k += m + 1;
	      break;
	    }
	}
	} //fin else
      nt++;
      if (nt >= *nuse) break;
    } //fin for
  } //fin for
} //fin

void Xtranslate(double *x, 
		int mdim, 
		int nrnodes, 
		int nsample, 
		int *bestvar, 
		int *bestsplit, 
		int *bestsplitnext,
		double *xbestsplit, 
		int *nodestatus, 
		int treeSize) {
/*
 this subroutine takes the splits on numerical variables and translates them
 back into x-values.  It also unpacks each categorical split into a 
 32-dimensional vector with components of zero or one--a one indicates that 
 the corresponding category goes left in the split.
*/

    int i, m;

    for (i = 0; i < treeSize; ++i) {
      if (nodestatus[i] == 1) {
	m = bestvar[i] - 1;
	xbestsplit[i] = 0.5 * (x[m + (bestsplit[i] - 1) * mdim] +
			       x[m + (bestsplitnext[i] - 1) * mdim]);
	
      }
    }
}

int pack(int nBits, int *bits) {
    int i = nBits, pack=0;
    while (--i >= 0) pack += bits[i] << i;
    return(pack);
}

void unpack(int pack, int *bits) {
/* pack is a 4-byte integer.  The sub. returns icat, an integer array of 
   zeroes and ones corresponding to the coefficients in the binary expansion 
   of pack. */   
    int i;
    for (i = 0; pack != 0; pack >>= 1, ++i) bits[i] = pack & 1;
}

#ifdef OLD

double oldpack(int l, int *icat) {
    /* icat is a binary integer with ones for categories going left 
     * and zeroes for those going right.  The sub returns npack- the integer */
    int k;
    double pack = 0.0;

    for (k = 0; k < l; ++k) {
	if (icat[k]) pack += R_pow_di(2.0, k);
    }
    return(pack);
}


void oldunpack(int l, int npack, int *icat) {
/*      
 * npack is a long integer.  The sub. returns icat, an integer of zeroes and
 * ones corresponding to the coefficients in the binary expansion of npack.
 */   
    int i;
    zeroInt(icat, 32);
    icat[0] = npack % 2; 
    for (i = 1; i < l; ++i) {
	npack = (npack - icat[i-1]) / 2;
	icat[i] = npack % 2;
    }
}



#endif /* OLD */
