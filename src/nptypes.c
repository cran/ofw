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
#include <R_ext/Utils.h>
#include "nptypes.h"

p_matrixdouble_t matrixdoubleMatrix(int nrows, int ncols) {

  int i;

  p_matrixdouble_t pmat = (p_matrixdouble_t) S_alloc(1, sizeof(matrixdouble_t));
  pmat->storage = (double *) S_alloc(nrows*ncols, sizeof(double));
  pmat->array = (double **) S_alloc(nrows, sizeof(double *));

  pmat->nrows = nrows;
  pmat->ncols = ncols; 
  
  for (i = 0; i < nrows; i++)
    {
      pmat->array[i] = pmat->storage + (i * ncols);
    }
  
  return(pmat);
}

p_matrixdouble_t matrixdoubleMatrixfromV(int nrows, int ncols, double* x){

  int i;

  p_matrixdouble_t pmat = (p_matrixdouble_t) S_alloc(1, sizeof(matrixdouble_t));

  pmat->storage = x;

  pmat->array = (double **) S_alloc(nrows, sizeof(double *));

  pmat->nrows = nrows;
  pmat->ncols = ncols; 
  
  for (i = 0; i < nrows; i++)
    {
      pmat->array[i] = pmat->storage + (i * ncols);
    }
  
  return(pmat);

}

void matrixdoublePrint(p_matrixdouble_t pmat) {
  
  int i,j;

  for (i = 0; i < pmat->nrows; i++)
    {
      for (j = 0; j < pmat->ncols; j++)
	printf("%5f ", pmat->array[i][j]);
      printf("\n");
    }
}

p_matrixint_t matrixMatrix(int nrows, int ncols) {

  int i;

  p_matrixint_t pmat = (p_matrixint_t) S_alloc(1, sizeof(matrixint_t));

  pmat->storage = (int *) S_alloc(nrows*ncols, sizeof(int));

  pmat->array = (int **) S_alloc(nrows, sizeof(int *));

  pmat->nrows = nrows;
  pmat->ncols = ncols; 
  
  for (i = 0; i < nrows; i++)
    {
      pmat->array[i] = pmat->storage + (i * ncols);
    }
  
  return(pmat);
}

p_matrixint_t matrixMatrixfromV(int nrows, int ncols, int* x) {

  int i;

  p_matrixint_t pmat = (p_matrixint_t) S_alloc(1, sizeof(matrixint_t));
  pmat->storage = x;
  pmat->array = (int **) S_alloc(nrows, sizeof(int *));

  pmat->nrows = nrows;
  pmat->ncols = ncols; 
  
  for (i = 0; i < nrows; i++)
    {
      pmat->array[i] = pmat->storage + (i * ncols);
    }
  
  return(pmat);
}

void matrixPrint(p_matrixint_t pmat) {
  
  int i,j;

  for (i = 0; i < pmat->nrows; i++)
    {
      for (j = 0; j < pmat->ncols; j++)
	printf("%5d ", pmat->array[i][j]);
      printf("\n");
    }
}

void matrixCopy(p_matrixint_t pmattarget, p_matrixint_t pmatsource) {
  memcpy(pmattarget->storage, pmatsource->storage, 
	 sizeof(int) * pmatsource->ncols * pmatsource->nrows);
}

