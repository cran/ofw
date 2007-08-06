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

#ifndef NPT_H
#define NPT_H

/** Matrix of double 
 *  
 */
typedef struct {
  int nrows; /**< Le nombre de lignes*/
  int ncols; /**< Le nombre de colonnes*/
  double* storage; /**< La zone de stockage*/
  double** array; /**< La zone d'indexation*/
} matrixdouble_t, *p_matrixdouble_t;

/** \brief The "constructor".
 * 
 * @param nrows Le nombre de lignes.
 * @param ncols Le nombre de colonnes.
 * @return Un pointeur sur une matrice de doubles.
 */
p_matrixdouble_t matrixdoubleMatrix(int nrows, int ncols);

/** \brief another "constructor".
 * 
 * @param nrows Le nombre de lignes.
 * @param ncols Le nombre de colonnes.
 * @param x Vector of double, the values
 * @return Un pointeur sur une matrice de doubles.
 */
p_matrixdouble_t matrixdoubleMatrixfromV(int nrows, int ncols, double* x);

/** \brief a raw print of the matrix type
 *
 * @param pmat the pointer on the matrix to print. 
 */
void matrixdoublePrint(p_matrixdouble_t pmat);

/** Type matrice d'entier. 
 *  
 */
typedef struct {
  int nrows; /**< Le nombre de lignes*/
  int ncols; /**< Le nombre de colonnes*/
  int* storage; /**< La zone de stockage*/
  int** array; /**< La zone d'indexation*/
} matrixint_t, *p_matrixint_t;

/** \brief Le "constructeur" du type matrice d'entiers.
 * 
 * @param nrows Le nombre de lignes.
 * @param ncols Le nombre de colonnes.
 * @return Un pointeur sur une matrice d'entiers.
 */
p_matrixint_t matrixMatrix(int nrows, int ncols);

/** \brief another "constructor".
 * 
 * @param nrows Le nombre de lignes.
 * @param ncols Le nombre de colonnes.
 * @param x Vector of int, the values
 * @return Un pointeur sur une matrice de doubles.
 */
p_matrixint_t matrixMatrixfromV(int nrows, int ncols, int* x);

/** \brief a raw print of the matrix type
 *
 * @param pmat the pointer on the matrix to print. 
 */
void matrixPrint(p_matrixint_t pmat);

/** \brief copy of matrix of int
 * Assuming they have the same size and that they are 
 * previously constructed.
 * @param pmattarget the pointer on the matrix.  
 * @param pmatsource the pointer on the matrix. 
 */
void matrixCopy(p_matrixint_t pmattarget, p_matrixint_t pmatsource);

#endif
