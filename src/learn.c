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
#include "learn.h"
#include "agregTree.h"
#include "classTree.h"



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
		) {
  

    int nsample0, mdim, nclass, mtry, ntest, nsample, ndsize, mimp, nimp, nuse, nrightall, keepInbag;
    int jb, j, n, m, k,  replace, stratify, trace, *nright, Ntree, i;

    int *out, *bestsplitnext, *bestsplit, *nodepop, *jin, *nodex,*nodexts, *nodestart, *ta, *ncase, *jerr, *varUsed, *jtr, *classFreq, *idmove, *jvr, *nind, *jts, *oobpair;


    p_matrixint_t p_m_at, p_m_a, p_m_b;

    p_matrixdouble_t p_m_x; /* x wrapped */

    p_matrixint_t p_m_counttr; /* counttr wrapped */
    p_matrixdouble_t p_m_countts; /* countts wrapped */

    int  last, ktmp, anyEmpty, ntry;

    double av=0.0;
    
    double *tgini, *tx, *wl, *classpop, *tclasspop, *win, *tp, *wr;

    int *index, Ngene, *omega, Weight;
    double *probaTemp, *erroob;
     
    mdim     = dimx[0];
    nsample0 = dimx[1];
    nclass   = (*ncl==1) ? 2 : *ncl;
    Ntree    = *ntree;   
    mtry     = *nvar;
    Ngene    = *ngene;
    ntest    = *nts; 
    Weight= *weight;
    nsample = nsample0;  


    tgini =      (double *) S_alloc(mdim, sizeof(double));
    wl =         (double *) S_alloc(nclass, sizeof(double));
    wr =         (double *) S_alloc(nclass, sizeof(double));
    classpop =   (double *) S_alloc(nclass* *nrnodes, sizeof(double));
    tclasspop =  (double *) S_alloc(nclass, sizeof(double));
    tx =         (double *) S_alloc(nsample, sizeof(double));
    win =        (double *) S_alloc(nsample, sizeof(double));
    tp =         (double *) S_alloc(nsample, sizeof(double));

    out =           (int *) S_alloc(nsample, sizeof(int));
    bestsplitnext = (int *) S_alloc(*nrnodes, sizeof(int));
    bestsplit =     (int *) S_alloc(*nrnodes, sizeof(int));
    nodepop =       (int *) S_alloc(*nrnodes, sizeof(int));
    nodestart =     (int *) S_alloc(*nrnodes, sizeof(int));
    jin =           (int *) S_alloc(nsample, sizeof(int));
    nodex =         (int *) S_alloc(nsample, sizeof(int));
    nodexts =       (int *) S_alloc(ntest, sizeof(int));
    ta =            (int *) S_alloc(nsample, sizeof(int));
    ncase =         (int *) S_alloc(nsample, sizeof(int));
    jerr =          (int *) S_alloc(nsample, sizeof(int));
    varUsed =       (int *) S_alloc(mdim, sizeof(int)); 
    jtr =           (int *) S_alloc(nsample, sizeof(int));
    jvr =           (int *) S_alloc(nsample, sizeof(int));
    classFreq =     (int *) S_alloc(nclass, sizeof(int));
    jts =           (int *) S_alloc(ntest, sizeof(int));
    idmove =        (int *) S_alloc(nsample, sizeof(int));

    p_m_x =  matrixdoubleMatrixfromV(nsample, mdim, x);

    p_m_countts =  matrixdoubleMatrixfromV(ntest, nclass, countts);

    p_m_counttr =  matrixMatrixfromV(nsample, nclass, counttr);

    p_m_at = matrixMatrix(nsample, mdim);
    p_m_a = matrixMatrix(nsample, mdim);
    p_m_b = matrixMatrix(nsample, mdim);
  
    nright =        (int *) S_alloc(nclass, sizeof(int));

    /* allocation */


     probaTemp =     (double *) S_alloc(mdim, sizeof(double));
     index =     (int *) S_alloc(mdim, sizeof(int)); 
     erroob=     (double *) S_alloc(Ntree, sizeof(double));
     omega =     (int *) S_alloc(mtry, sizeof(int)); 


     //initialisation des variables pr omega
     zeroInt(index, mdim);
     zeroDouble(probaTemp, mdim);
     
     //sort the P-weighted variables
     for (i=0; i<mdim;i++) {
       index[i]=i;
       probaTemp[i]=proba[i];
     }
     revsort(probaTemp, index, mdim);
     
     /*
     for (i=0; i<Ngene; i++){
       printf("index  %d ", index[i] );
     }
     */
     /*for (i=0; i<nclass; i++){
       printf("classweight  %lf ", classWeight[i] );
       }*/


     //BIG LOOP ON TREES

    /*    INITIALIZE FOR RUN */
     zeroDouble(p_m_countts->storage, ntest * nclass);
     zeroInt(p_m_counttr->storage, nclass * nsample);
     
     zeroDouble(tgini, mdim);
     zeroDouble(erroob,Ntree);  
     zeroDouble(errin,Ntree); 
     zeroDouble(errts,Ntree); 

     zeroInt(outcl,nsample); 
     zeroInt(outclts,ntest); 
     
     makeA(p_m_x ,
	   p_m_at, 
	   p_m_b); 
     
     R_CheckUserInterrupt();
     GetRNGstate(); 
     

     zeroInt(ndbigtree, Ntree); 
     
     //boucle sur les arbres
     for(jb = 0; jb < Ntree; jb++) {
      do {

	zeroInt(nodestatus , *nrnodes);
	zeroInt(treemap, 2 * *nrnodes);
	zeroDouble(xbestsplit, *nrnodes);
	zeroInt(nodeclass, *nrnodes);
	zeroInt(varUsed, mdim);
	//initialisations en plus
	zeroInt(bestvar, *nrnodes);
	zeroInt(ncase, nsample);
		
	zeroDouble(tclasspop, nclass);  //je reinitialise
	zeroInt(jin, nsample);
	zeroDouble(win, nsample);
	zeroInt(jts, ntest);
	zeroInt(nodexts, ntest);
	zeroInt(jtr, nsample);
	zeroInt(nodex, nsample);

	
	for (n = 0; n < nsample; ++n) { 
	  tclasspop[cl[n] - 1] += classwt[cl[n]-1];
	  win[n] += classwt[cl[n]-1];
	  jin[n] = 1;
	}
	  	
	//trouver le min pour nodesize
	double minSize;
	minSize= tclasspop[nclass-1];
	for (n=0; n<(nclass-1); ++n){
	  minSize=(tclasspop[n]<=minSize) ? tclasspop[n]:minSize;
	}
	ndsize=(int) (minSize);  
	
	
	/* Copy the original a matrix back. */
	
	matrixCopy(p_m_a, p_m_at);
	
	modA(p_m_a, 
	     &nuse,
	     ncase, 
	     jin);
	
	
	/* draw the mtry variables among the Ngene variables (in omega)*/   
	zeroInt(omega, mtry);   //mtry=nvar is the variables number in omega 
	if (Ngene <=mtry){
	  for (i=0; i<mtry; i++){
	    omega[i]=index[i];
	  }
	}
	else{
	  for (i=0; i<mtry; i++){
	    k = unif_rand() * Ngene;
	    omega[i]=index[k];
	  }
	}



	/*
	for (i=0; i<mtry; i++){
	   printf("omega  %d ", omega[i] );
	   }*/



	/* build the tree with the omega variables */
	buildtree(p_m_a->storage, 
		  p_m_b->storage, 
		  cl,  
		  &mdim, 
		  &nsample, 
		  &nclass, 
		  treemap, 
		  bestvar,
		  bestsplit, 
		  bestsplitnext, 
		  tgini, 
		  nodestatus, 
		  nodepop, 
		  nodestart, 
		  classpop, 
		  tclasspop,
		  ta, 
		  nrnodes, 
		  idmove, 
		  &ndsize, 
		  ncase,   
		  &mtry, 
		  varUsed, 
		  nodeclass, 
		  ndbigtree + jb, 
		  win, 
		  wr, 
		  wl,
		  &nuse,
		  omega
		  );
	/* if the "tree" has only the root node, start over */
	
      } while (ndbigtree[jb] == 1);
      
      Xtranslate(p_m_x->storage, mdim, *nrnodes, nsample, bestvar, 
		 bestsplit, bestsplitnext, xbestsplit,
		 nodestatus,
		 ndbigtree[jb]);
      
      /*  Get test set error */
	
		predictClassTree(xts, 
			 ntest, 
			 mdim, 
			 treemap,
			 nodestatus, 
			 xbestsplit,
			 bestvar, 
			 nodeclass,
			 jts, 
			 nodexts
			 );
			 
	
	TestSetError(p_m_countts, jts, clts, outclts, ntest, nclass, jb+1,
		     errts + jb , cut, Weight, classWeight);
		     
     
      //printf("errts %lf", errts[jb]);
      
      /*  Get out-of-bag predictions and errors. */

      predictClassTree(p_m_x->storage, 
		       nsample, 
		       mdim, 
		       treemap,
		       nodestatus, 
		       xbestsplit,
		       bestvar, 
		       nodeclass,
		       jtr, 
		       nodex
		       );
      

      /* Compute the inbag  error rate. */
      inbagError(nsample, nclass, cl, p_m_counttr, errin +jb, outcl, cut, jtr, jb+1, Weight, classWeight);
      
      //printf("errIn %lf \n ", errin[jb]);
      
      /*
      if(jb==(Ntree-1)) {
	printf("errts %lf  ", errts[jb]);
	printf("errIn  %lf \n", errin[jb]);

	for (n=0; n<ntest; n++) printf("jts  %d \n", jts[n]);
	for (n=0; n<ntest; n++) printf("outclts  %d \n", outclts[n]);
	
	}*/ 

      
      
    } //fin jb
     

	
	R_CheckUserInterrupt();
#ifdef win32
	R_ProcessEvents();
#endif
	
	PutRNGstate();

}  //fin classError





///////////////////////////////////////////////////////
/*
  Modified by A. Liaw 1/10/2003 (Deal with cutoff)
  Re-written in C by A. Liaw 3/08/2004
  puis Kim-Anh
*/
void inbagError(
	 int nsample, 
	 int nclass,
	 int *cl, 
	 p_matrixint_t p_m_counttr,    //int *counttr,  
	 double *errtr, //taux d erreur (sortie)
	 int *jest,  //equiv outcl: classe predite
	 double *cutoff,
	 int *jtr,
	 int nvote,
	 int weight,
	 double *classWeight
	 ) {
  // calcul le taux d erreur oob et renvoie errtr
  int j, n;
  double qq, smax, smaxtr, sumWeight;
  sumWeight=0.0;
  
  //increment the predicted votes
  for (n = 0; n < nsample; ++n) p_m_counttr->array[n][jtr[n] - 1] ++; 
      
  for (n = 0; n < nsample; ++n) {
    smax = 0.0;
    smaxtr = 0.0;
    for (j = 0; j < nclass; ++j) {
      qq = (((double) p_m_counttr->array[n][j]) / nvote) / cutoff[j]; 
      if (j+1 != cl[n]) smax = (qq > smax) ? qq : smax;
      /* if vote / cutoff is larger than current max, re-set max and 
	 change predicted class to the current class */
      if (qq > smaxtr) {
	smaxtr = qq;
	jest[n] = j+1;  //num de la classe de l'obs
      }
      /* break tie at random */
      if (qq == smaxtr && unif_rand() > 0.5) {
	smaxtr = qq;
	jest[n] = j+1;
      }
    }
    if (jest[n] != cl[n]) {
      if (weight==0) {errtr[0] += 1.0;}
     else {errtr[0] += 1.0 * classWeight[cl[n]-1];}
    }
    if (weight==1) {sumWeight+=classWeight[cl[n]-1];}
  } //fin for sur n
  if (weight==0){errtr[0] /= nsample;}
  else { errtr[0] /= sumWeight;}

}




void TestSetError(
		  p_matrixdouble_t p_m_countts, 
		  int *jts, 
		  int *clts, 
		  int *jet, 
		  int ntest, 
		  int nclass, 
		  int nvote, 
		  double *errts,  
		  double *cutoff,
		  int weight,
		  double *classWeight
		  ) {

  int j, n;
  double cmax, crit, sumWeight;
  sumWeight=0.0;

   for (n = 0; n < ntest; ++n) p_m_countts->array[n][jts[n]-1] += 1.0;

  //  Prediction is the class with the maximum votes 
  for (n = 0; n < ntest; ++n) {
    cmax=0.0;
    for (j = 0; j < nclass; ++j) {
      crit = (p_m_countts->array[n][j] / nvote) / cutoff[j];
      if (crit > cmax) {
	jet[n] = j+1;
	cmax = crit;
      }
      //  Break ties at random: 
      if (crit == cmax && unif_rand() > 0.5) {
		jet[n] = j+1;
		cmax = crit;
	    }
	}
  }

  for (n = 0; n < ntest; ++n) {
    if (jet[n] != clts[n]) {
      if (weight==0) {errts[0] += 1.0;}
      else {errts[0] += 1.0 * classWeight[clts[n]-1]; }
    }
    if (weight==1) {sumWeight+=classWeight[clts[n]-1];}
    //printf("sumWeight  %lf \n", sumWeight);} 

  }
  if (weight==0){errts[0] /= ntest;}
  else { errts[0] /= sumWeight;}
  //printf("sumWeight  %lf \n", sumWeight);
  
  //errts[0] /= ntest;
       
}

