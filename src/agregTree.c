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
#include "agregTree.h"
#include "classTree.h"

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
	     ) {
    /******************************************************************
     *  C routine for Optimal Feature Weighting algorithù:  get input from R and drive the C routines.
     *
     *  Input:
     *  x:        matrix of predictors (transposed!)
     *  dimx:     two integers: number of variables and number of cases
     *  cl:       class labels of the data
     *  ncl:      number of classes in the response
     * Options:   4 integers: (0=no, 1=yes)    
     *     trace:     how often to check the stopping criterion? 
     *     keepInbag: keep the last forest check the construction of the trees
     *     keepf:     keep the last forest
     *     weight:    weight the classes when calculating the internal error
     *
     *  ntree:    number of trees
     *  nforest:  number of iterations of the algorithm (= number of forests)
     *  nvar:     number of predictors to use for each tree
     * nstable:   number of stable weighted variables expected
     * nrnodes:   maximum number of nodes in each tree
     * ndbigtree  
     * bestvar
     *treemap
     *nodeclass: class of the node of each tree
     *xbestsplit
     *
     *Output:
     * maxiter:   number of iterations reached
     *inbag:      observations used to construct each tree for the last iteration
     *moyenne:    internal error rate for each iteration
     *proba:      probability weight on each variable
     ******************************************************************/

    int mdim, nclass, mtry, ntest, nsample, ndsize, nuse, keepInbag, keepf;
    int jb, j, n, m, k, idxByNnode, idxByNsample, trace, Ntree, Nforest, weight, Nstable, Maxiter;

    int *out, *bestsplitnext, *bestsplit, *nodepop, *jin, *nodex, *nodestart, *ta, *ncase, *jerr, *varUsed,*jtr, *classFreq,  *idmove;

    double *probaTemp, *classwt;
    int *index1, *index2;
    double *gradient, *erroob;
    int *omega;
    int i;
    double *sumProba;
    int *perm;

    p_matrixint_t p_m_at, p_m_a, p_m_b;

    p_matrixdouble_t p_m_x; /* x wrapped */

    int anyEmpty;

    
    double *tgini, *tx, *wl, *classpop, *tclasspop, *win, *tp, *wr;

    trace    = Options[0];                                   
    keepInbag = Options[1]; 
    keepf = Options[2];  
    weight= Options[3]; 
    mdim     = dimx[0];
    nsample = dimx[1];
    nclass   = (*ncl==1) ? 2 : *ncl;
    Ntree    = *ntree;  
    Nforest = *nforest;
    mtry     = *nvar;
    Nstable = *nstable;
                   
    if (trace == 0) trace = Nforest+1;

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
    ta =            (int *) S_alloc(nsample, sizeof(int));
    ncase =         (int *) S_alloc(nsample, sizeof(int));
    jerr =          (int *) S_alloc(nsample, sizeof(int));
    varUsed =       (int *) S_alloc(mdim, sizeof(int)); 
    jtr =           (int *) S_alloc(nsample, sizeof(int));
    classFreq =     (int *) S_alloc(nclass, sizeof(int));
    classwt =     (double *) S_alloc(nclass, sizeof(double));
    idmove =        (int *) S_alloc(nsample, sizeof(int));

    p_m_x =  matrixdoubleMatrixfromV(nsample, mdim, x);

    p_m_at = matrixMatrix(nsample, mdim);
    p_m_a = matrixMatrix(nsample, mdim);
    p_m_b = matrixMatrix(nsample, mdim);
  
    /* allocation */
     omega =        (int *) S_alloc(mtry, sizeof(int));
     gradient =     (double *) S_alloc(mdim, sizeof(double));

     probaTemp =     (double *) S_alloc(mdim, sizeof(double));
     index1 =     (int *) S_alloc(mdim, sizeof(int)); 
     index2 =     (int *) S_alloc(mdim, sizeof(int));

     erroob=     (double *) S_alloc(Ntree, sizeof(double));
     sumProba =     (double *) S_alloc(mdim, sizeof(double));
     perm =     (int *) S_alloc(mdim, sizeof(int)); 

     /*for the class weights*/
     int *numWeight, *inside;
     double *classWeight; 
     numWeight =     (int *) S_alloc(nclass, sizeof(int)); 
     classWeight =     (double *) S_alloc(nclass, sizeof(double)); 
     inside =     (int *) S_alloc(nsample, sizeof(inside)); 
    
     /*needed for the stopping criterion */
     zeroInt(index1, mdim);
     zeroInt(index2, mdim);
     zeroDouble(probaTemp, mdim);
     
     /* initialize with a uniform probability */
    proba_init(proba, &mdim);
    
    /* initialisation of gradient and internal error rate */
    zeroDouble(gradient, mdim);
    zeroDouble(moyenne, Nforest);
    
    //loop for Nforests
    for (n = 1; n <=Nforest ; n++)
      {
	if(n%trace ==0) printf(" iteration \%5d ",n);
	if ((trace ==0) && (n==Nforest))  printf(" iteration \%5d ",n);
	
	classForest(p_m_x, 
		    dimx, 
		    cl, 
		    ncl, 
		    Options, 
		    ntree, 
		    nforest,
		    nvar,
		    classwt,  
		    nrnodes, 
		    ndbigtree, 
		    nodestatus, 
		    bestvar, 
		    treemap, 
		    nodeclass, 
		    xbestsplit, 
		    erroob, 
		    inbag,
		    tgini,
		    wl,
		    wr,
		    classpop,   
		    tclasspop, 
		    tx, 
		    win,
		    tp,
		    out,
		    bestsplitnext,
		    bestsplit,
		    nodepop,
		    nodestart,
		    jin,
		    nodex,
		    ta,
		    ncase,         
		    jerr,
		    varUsed,
		    jtr,
		    classFreq ,
		    idmove,     
		    p_m_at,
		    p_m_a,
		    p_m_b,
		    proba,
		    omega,
		    gradient,
		    n,
		    moyenne,
		    sumProba,
		    perm,
		    inside,
		    numWeight,
		    classWeight
		    );
	/*gradient descent */
	for (i =0; i<mdim; i++){
	  proba[i]=proba[i] - gradient[i]/(1+n);
	}
	
	/*Project the probability on the simplex */
	simplexe(proba, &mdim);
	

     	/*A criterion to stop the forest when the weighted variables are the same */	
	if(n==0){
	  zeroDouble(probaTemp, mdim);  
	  for (i=0; i<mdim;i++) {
	    index1[i]=i;
	    probaTemp[i]=proba[i];
	  }
	  revsort(probaTemp, index1, mdim);
	}
	
	if (n%trace==0) {
	  zeroDouble(probaTemp, mdim); 
	  zeroInt(index2, mdim);
	  for (i=0; i<mdim;i++) {
	    index2[i]=i;
	    probaTemp[i]=proba[i];
	  }
	  revsort(probaTemp, index2, mdim);
	  
	  /* compare if the first weighted variables are the same every trace */
	  int K;
	  K=0;
	  for(i=0; i<Nstable; i++) {
	  for(j=0; j<Nstable; j++){
	    if(index1[i]==index2[j]) K++;
	  }
	}
	  
	  if (K==Nstable) {
	    *maxiter=n;
	    n=Nforest;
	  }
	  printf(" stable variables: \%5d ",K); 
	  printf(" \n");
	  
	  for (i=0; i<mdim; i++) index1[i]=index2[i];
	}
	
		
      }  //fin Nforest
      printf("\n");
      //printf("iteration \%5d \n",n-1);
}  //fin classAgregTree







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
		 double *erroob, 
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
		 int * omega,
		 double * gradient,
		 int iter,
		 double * moyenne,
		 double * sumProba,
		 int * perm,
		 int * inside,
		 int * numWeight,
		 double * classWeight
		 ) {
 
  int  mdim, nclass, mtry, ntest, nsample, ndsize, nuse, keepInbag, keepf;
  int jb, j, n, m, k, idxByNnode, idxByNsample, trace, Ntree, Nforest, weight;
  int i;

  int  anyEmpty;
   
  
  trace    = Options[0];                       
  keepInbag = Options[1];
  keepf= Options[2];
  weight= Options[3];
  mdim     = dimx[0];
  nsample = dimx[1];
  nclass   = (*ncl==1) ? 2 : *ncl;
  Ntree    = *ntree;
  Nforest = *nforest;
  mtry     = *nvar;
  

    /* Count number of cases in each class. */
    zeroInt(classFreq, nclass);
    zeroDouble(classwt, nclass);
    for (n = 0; n < nsample; ++n) classFreq[cl[n] - 1] ++;
    /* Normalize class weights. */
    for (i = 0; i < nclass; ++i) {
      classwt[i] = ((double) classFreq[i]) / nsample;
     }
     
     zeroInt(numWeight, nclass);
     zeroDouble(classWeight, nclass);
     if (weight==1){
	  //compute the weights
	  for (i =0; i<nclass; i++){
	    classWeight[i]=(double)1.00/(classFreq[i] * nclass);
	  }
	}//fin weight
     
  

    /*    INITIALIZE FOR RUN */
    
    zeroDouble(tgini, mdim);
    zeroDouble(erroob,Ntree); 

    makeA(p_m_x ,
	  p_m_at, 
	  p_m_b); 

    R_CheckUserInterrupt();

    /* Starting the main loop over number of trees. */
    GetRNGstate(); 

    idxByNnode = 0;
    idxByNsample = 0;
    zeroInt(ndbigtree, Ntree); 


    //loop on the trees
    for(jb = 0; jb < Ntree; jb++) {
      do {
	zeroInt(nodestatus + idxByNnode, *nrnodes);
	zeroInt(treemap + 2*idxByNnode, 2 * *nrnodes);
	zeroDouble(xbestsplit + idxByNnode, *nrnodes);
	zeroInt(nodeclass + idxByNnode, *nrnodes);
	zeroInt(varUsed, mdim);
	zeroInt(bestvar + idxByNnode, *nrnodes);
	zeroInt(ncase, nsample);

	//zeroInt(numWeight, nclass);
	//zeroDouble(classWeight, nclass);
	   //zeroDouble(sampleWeight, nsample);	
	
	//draw a bootstrap sample with no empty class
	anyEmpty=1;
	int oob;
	do {
	  anyEmpty = 0;
	  oob=0;
	  zeroDouble(tclasspop, nclass); 
	  zeroInt(jin, nsample);
	  zeroDouble(win, nsample);
	  zeroInt(inside, nsample);
	  
	  for (n = 0; n < nsample; ++n) {  
	    k = unif_rand() * nsample;
	    inside[k] +=1; //inside counts how many times a sample is drawn
	    tclasspop[cl[k] - 1] += classwt[cl[k]-1];
	    win[k] += classwt[cl[k]-1];
	    jin[k] = 1;
	  }
	  
	  
	  /* check if any class is missing in the sample */
	  for (n = 0; n < nclass; ++n) {
	    if (tclasspop[n] == 0) anyEmpty = 1;
	    if (jin[n] == 0) oob++;
	  }
	} while (anyEmpty || oob==0); 


/*	//compute the weights
	//int sum=0.0;
	if (weight==1){
	  //how many samples in class i ?
	  for(i =1; i<=nclass; i++){
	    for (n=0; n<nsample; n++){
	      if(cl[n]==i) numWeight[i-1]+= inside[n];  //classes are indexed from 1
	      
	    }
	  }
	  //compute the weights
	  for (i =0; i<nclass; i++){
	     classWeight[i]=(double)1.00/(numWeight[i] * nclass);
	    //classWeight[i]=(double)nsample/numWeight[i];
	    //sum +=classWeight[i];
	    printf("classWeight  %f \n", classWeight[i]);
	    printf("numWeight  %f \n", numWeight[i]);
	  }
	  
	  //for (i =0; i<nclass; i++) classWeight[i]/=sum;

	}//fin weight   */
	
	
	
	//determine the minimal size of each node
	double minSize;
	minSize= tclasspop[nclass-1];
	for (n=0; n<(nclass-1); ++n){
	  minSize=(tclasspop[n]<=minSize) ? tclasspop[n]:minSize;
	}
	ndsize=(int) (minSize);  

	
	/* If need to keep indices of inbag data, do that here. */
	if (keepInbag) {
	  for (n = 0; n < nsample; ++n) {
	    inbag[n + idxByNsample] = inside[n];  //au lieu de jin[n]
	  }
	}
	
	/* Copy the original a matrix back. */

	matrixCopy(p_m_a, p_m_at);
	
	modA(p_m_a, 
	     &nuse,
	     ncase, 
	     jin);
	
	
	/* draw the omega variables */

	zeroInt(omega, mtry);
	/*draw  mtry omega variables with probability proba */
	draw_omega(proba,&mdim,omega, &mtry, sumProba,perm); 

	/* build the tree with the omega variables */
	buildtree(p_m_a->storage, 
		  p_m_b->storage, 
		  cl,  
		  &mdim, 
		  &nsample, 
		  &nclass, 
		  treemap + 2*idxByNnode, 
		  bestvar + idxByNnode,
		  bestsplit, 
		  bestsplitnext, 
		  tgini, 
		  nodestatus + idxByNnode, 
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
		  nodeclass + idxByNnode, 
		  ndbigtree + jb, 
		  win, 
		  wr, 
		  wl,
		  &nuse,
		  omega
		  );
	/* if the "tree" has only the root node, start over */
      } while (ndbigtree[jb] == 1);
      
      Xtranslate(p_m_x->storage, mdim, *nrnodes, nsample, bestvar + idxByNnode, 
		 bestsplit, bestsplitnext, xbestsplit + idxByNnode,
		 nodestatus + idxByNnode,
		 ndbigtree[jb]);
      
      
      /*  Get out-of-bag predictions and errors. */

      predictClassTree(p_m_x->storage, 
		       nsample, 
		       mdim, 
		       treemap + 2*idxByNnode,
		       nodestatus + idxByNnode, 
		       xbestsplit + idxByNnode,
		       bestvar + idxByNnode, 
		       nodeclass + idxByNnode,
		       jtr, 
		       nodex
		       );


      //compute the oob error rate (erroob) on each tree
      if (weight==1){oobTreeWeight(nsample, cl, jtr, jin, erroob, jb, classWeight);
      }else{
	oobTree(nsample, cl, jtr, jin, erroob, jb);
      }
      // printf("ooberror  %lf \n", erroob[jb] );

      //compute the gradient for each tree and each drawn omega 
      for (i=0; i<mtry; i++){
	gradient[omega[i]] += erroob[jb]/(10000*proba[omega[i]]*(iter+1)) ;	
      }

      R_CheckUserInterrupt();
#ifdef win32
      R_ProcessEvents();
#endif
      if (keepf) idxByNnode += *nrnodes;
      if (keepInbag) idxByNsample += nsample;


    } //fin jb


    /* get the mean error of all trees in the forest */
    for (jb=0; jb<Ntree; jb++){
      moyenne[iter] +=  erroob[jb];
    } 
    moyenne[iter] /= Ntree; 
    //printf("error  %lf \n",moyenne[iter]);
    
    PutRNGstate();
    

}  //fin classAgregTree




/////////////////////////////////////////////////////////////////
void oobTree(
	 int nsample, 
	 //int nclass, 
	 int *cl, 
	 int *jtr,  
	 int *jin, 
	 double *errtr, //taux d erreur (sortie)
	 int jb
	 ) {
  // calcul le taux d erreur oob sur un arbre et renvoie errtr
  int  n, noob;
  
  noob = 0;
  for (n = 0; n < nsample; ++n) {
    if (jin[n]==0) {  //si lobs est oob
      noob++;
      if (jtr[n] != cl[n]) errtr[jb] += 1.0;
    }
  }

  //if( noob==0){
  //printf("noob %d \n", noob);
  errtr[jb] /= noob;
  // printf("errtr %lf \n", errtr[jb]);
  // if( noob==0) printf("errtr %lf \n", errtr[jb]);
}



/////////////////////////////////////////////////////////////////
void oobTreeWeight(
	 int nsample, 
	 int *cl, 
	 int *jtr,  
	 int *jin, 
	 double *errtr, //error rate
	 int jb,
	 double *classWeight
	 ) {
  /* compute the oob error rate and returns errtr*/
  int  n;
  double noob;
  
  noob = 0.0;
  for (n = 0; n < nsample; ++n) {
    if (jin[n]==0) {  //if sample is oob
      noob += classWeight[cl[n]-1];
	if (jtr[n] != cl[n]) errtr[jb] += 1.0* classWeight[cl[n]-1];
    }
  } 
  
  errtr[jb] /= noob;
  //if( noob==0) printf("errtr %lf \n", errtr[jb]);
}


///////////////////////////////////////////////////////
//from Liaw
void predictClassTree(double *x,
		      int n, //number of tested samples
		      int mdim, 
		      int *treemap,
		      int *nodestatus, 
		      double *xbestsplit,
		      int *bestvar, 
		      int *nodeclass,
		      int *jts, 
		      int *nodex //node when the test sample lands
		      ) {
  //each tested sample runs down the tree, get the predicted class and the node number when it landed
    int m, i, k;
    //added: initialisation de jts et nodex
    zeroInt(jts, n);
    zeroInt(nodex, n);

    for (i = 0; i < n; ++i) {
	k = 0;
	while (nodestatus[k] != NODE_TERMINAL) {
	  m = bestvar[k] - 1;
	  /* Split by a numerical predictor */
	  // si valeur du tableau pr la variable m est inf au seuil
	  // on assigne a k le noeud courant (ncur) a gche sinon a droite
	  //ceci juska arriver au noeud final
	  k = (x[m + i * mdim] <= xbestsplit[k]) ?
	    treemap[k * 2] - 1 : treemap[1 + k * 2] - 1;
	 
	}
	/* Terminal node: assign class label */
	//get the class of k
	jts[i] = nodeclass[k];
	nodex[i] = k + 1;
    }
}




/* Initilize the uniform probability */
void proba_init(double *proba, int *mDim){ 
  int i;
  zeroDouble(proba, *mDim);
  for (i = 0; i < *mDim; ++i) proba[i]=1.0/ *mDim;
}


/* inspired from ProbSampleReplace from base package */
void draw_omega(double *proba,   //draw the omega variables with this probability
		int *mDim,        //number of total variables
		int *omega,       //subset of variables
		int *s_omega,      //subset size
		double *sumProba,
		int *perm
 		){   
  int i, j;
  double u;

  zeroDouble(sumProba, *mDim);
  zeroInt(perm, *mDim); 

  for (i =0; i<*mDim; i++) {
    sumProba[i]=proba[i];
    perm[i] = i;
  }
  
  // sort the probabilities into descending order 
  revsort(sumProba, perm, *mDim);

  // compute cumulative probabilities 
  for (i = 1 ; i < *mDim; i++)
    sumProba[i] += sumProba[i - 1];
  
  // compute the sample 
  for (i = 0; i < *s_omega; i++) {
    u = unif_rand();
    for (j = 0; j <(*mDim-1); j++) {
      if (u <= sumProba[j])
		break;
    }
    omega[i] = perm[j];
  }

  
} 


/*project the probability on the simplex */
void simplexe( double *proba,
	       int *mDim)
{
  int i,J,label;
  double somme,lambda;
  somme=0;

  for(i=0; i<*mDim ;i++){
    somme += proba[i];
  }
  lambda=(somme-1)/((double)*mDim);

  for(i=0;i< *mDim;i++){
    proba[i]=proba[i]-lambda;
  }

  label=0;
  while(label==0){
    label=1;
    somme=0;J=0;
    /* compute negative probabilities and set them to zero */
    for(i=0;i<*mDim;i++){
      if(proba[i]<=0){proba[i]=0;J++;}
      else{somme+=proba[i];}
    }
    /* project on the simplex the positive porbabilities */
    for(i=0;i<*mDim;i++){
      if(proba[i]>0){ proba[i]=proba[i]+(1-somme)/((double) (*mDim -J));}
      if(proba[i]<0){label=0;}
    }
  }  //fin while

}


