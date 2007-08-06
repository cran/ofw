/*******************************************************************
   Copyright (C) 2001-4 Leo Breiman, Adele Cutler and Merck & Co., Inc.
   Copyright (C) 2007 Kim-Anh Lê Cao, Patrick Chabrier, INRA,
   French National Institute for Agricultural Research.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.                            
*******************************************************************/

#include <R.h>
#include <Rmath.h>
#include "nputils.h"
#include "nptypes.h"

#include "classTree.h"


void buildtree(int *a, 
	       int *b, 
	       int *class, 
	       int *mDim, 
	       int *nsample,
	       int *nClass, 
	       int *treemap,  //
	       int *best,  //equiv de bestvar, !=bestVar !!
	       int *bestSplit, //
	       int *bestSplitNext,//
	       double *tgini, //
	       int *nodeStatus, //
	       int *nodePop, //
	       int *nodeStart,  
	       double *tclassCount,//equiv de classpop: tableau
	       double *classCount,//equiv de tclasspop (warning!):vecteur
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
	       //int *mIndex,
	       int *omega
	       ){
  /*
    attention : au nom des variables
    a SAVOIR: 
    -nrnodes est une valeur par defaut =2*(nsample+1) dans newpack.default.R
    -nuse: calculé ds modA : compte le nombre de inBag cad nb d'obs ds ech. bootstrap
    
    c     Buildtree consists of repeated calls to two subroutines, Findbestsplit
    c     and Movedata.  Findbestsplit does just that--it finds the best split of
    c     the current node.  Movedata moves the data in the split node right and
    c     left so that the data corresponding to each child node is contiguous.
    c     The buildtree bookkeeping is different from that in Friedman's original
    c     CART program.  ncur is the total number of nodes to date.
    c     nodestatus(k)=1 if the kth node has been split.  nodestatus(k)=2 if the
    c     node exists but has not yet been split, and =-1 of the node is terminal.
    c     A node is terminal if its size is below a threshold value, or if it is
    c     all one class, or if all the x-values are equal.  If the current node k
    c     is split, then its children are numbered ncur+1 (left), and
    c     ncur+2(right), ncur increases to ncur+2 and the next node to be split is
c     numbered k+1.  When no more nodes can be split, buildtree returns to the
c     main program.
  */
  
  int splitVar; 
  int ncur, j, kbuild, n, nc, k, kn, pp, m;
  int ndStart, ndEnd, splitStatus, ndendl, nbest;
  double decGini, popt1, popt2;
  double tieRand; // cas du noeud terminal on revient a la boucle kbuild pr autre split


  splitVar = 0;
  zeroInt( nodeStatus, *nrnodes);
  zeroInt( nodeStart, *nrnodes);
  zeroInt( nodePop, *nrnodes);
  zeroDouble( tclassCount, *nClass * *nrnodes);  

  //au noeud 1 on compte le nb d'obs par classe
  for (j= 0; j<*nClass; ++j) tclassCount[j]= classCount[j]; 

  //ncur : current node
  ncur = 1;

  nodeStart[0] = 1;  

  //nuse est le nb d obs presentes au noeud (cf fonction modA)
  //on a dc ds nodepop(0) le nb d'obs pr le noeud 1
  nodePop[0] = *nuse; 
  //le noeud 1 existe mais n est pas encore divise
  nodeStatus[0] = 2;   //rq idem

  
  //boucle sur kbuild
  for (kbuild=1; kbuild <= *nrnodes; ++kbuild){
    if ( kbuild > ncur) break; // goto 50
    if ( nodeStatus[kbuild - 1] != 2)  continue; //noeud splitte ou terminal

   //initialize for next call to findbestsplit
   //on definit le debut et la fin des endroits ou on peut faire un split
    ndStart = nodeStart[kbuild - 1 ];
    ndEnd = ndStart + nodePop[kbuild - 1] -1;
    //on compte le nb d obs par classe pr le noeud kbuild
    for (j = 0; j< *nClass; ++j){
      classCount[j]= tclassCount[j + (kbuild - 1) * *nClass];
    }
    splitStatus = 0;
    
    /*c RAPPEL: findbestsplit renvoie: 
      c             *(nbest): l'endroit ou splitter entre ndstart et ndend
      c             *splitVar (msplit): la variable ou splitter
      c             *splitStatus (jstat): le statut du noeud (qui peut etre -1 si noeud terminal)
      *decGini (decsplit) la decroissance de gini
      */
    findbestsplit(a,
		  b,
		  class, 
		  mDim,
		  nClass, 
		  &ndStart, 
		  &ndEnd, 
		  classCount, 
		  &splitVar, 
		  &decGini, 
		  &nbest, 
		  &splitStatus,
		  mtry, 
		  weight, 
		  wr,
		  wl,
		  omega
		  );   
    
    // cas du noeud terminal on revient a la boucle kbuild pr autre split:
    if( splitStatus == -1){
      nodeStatus[kbuild - 1] = -1;
      continue;   //revient au debut de la boucle
    } 
    else{
      //sinon on stock bestvar pr le noeud kbuild
      //rq splitVar = mvar+1 dc va de 1 à mdim
      best[kbuild-1]= splitVar; 
      varUsed[splitVar - 1] = 1; 
      //pas de decroissance negative, remise a zero
      if (decGini < 0.0) decGini= 0.0;
      //calcul de la decroissance de gini pr ce split
      tgini[splitVar - 1]+= decGini;
       
      //on stock l'obs splittee qui sera a gche et celle d apres qui sera a drte
      bestSplit[kbuild - 1]= a[(splitVar -1) + (nbest - 1) * *mDim];
      bestSplitNext[kbuild - 1]=a[(splitVar -1) + (nbest) * *mDim];
    } //fin if 
    /*RAPPEL: movedata reordonne a et renvoie ndendl=endroit ou splitter ainsi que ncase */
    
    movedata(a, ta, mDim, nsample, &ndStart, &ndEnd, idmove, ncase, &splitVar, &nbest, &ndendl);
	
    /*
      c     leftnode no.= ncur+1, rightnode no. = ncur+2.
      c     la division a lieu on separe les populations en 2
      c     ainsi que les indices de depart et arrivee pr chaque child node
    */
    nodePop[ncur]=ndendl - ndStart + 1;
    nodePop[ncur + 1]= ndEnd - ndendl;
    nodeStart[ncur]= ndStart;
    nodeStart[ncur + 1]= ndendl +1;
    
    /*find class populations in both nodes
      c     on recup a chaque noeud le nb d obs par classe (+ le poids)
      c        a gauche */
    for ( n = ndStart; n <= ndendl; ++n){
      nc = ncase[n-1]; 
      m=class[nc - 1];
      tclassCount[m-1 + (ncur)* *nClass]+=weight[nc-1];
    }
    //a droite
    for (n = ndendl+1; n<= ndEnd; ++n){
      nc = ncase[n-1];
      m=class[nc - 1];
      tclassCount[m-1 + (ncur + 1) * *nClass]+=weight[nc-1];
    }
    
    //check on nodestatus
    //on declare les nveaux child node comme non splitte
    nodeStatus[ncur]=2;
    nodeStatus[ncur + 1]=2;
    //si le nb d obs < ndsize on declare le noeud terminal
    if (nodePop[ncur] <= *ndsize) nodeStatus[ncur]= -1;
    if (nodePop[ncur+1] <= *ndsize) nodeStatus[ncur+1]= -1;
    
    //on compte le nb d'obs dans chaque noeud
    popt1=0.0;
    popt2=0.0;
    for (j =0; j < *nClass; ++j){
      popt1 += tclassCount[j + (ncur) * *nClass];
      popt2 += tclassCount[j + (ncur + 1) * *nClass];
    }
    
    //s'il n y a qu'1 classe ds un noeud on le declare comme terminal
    for (j =0; j < *nClass; ++j){
      if (tclassCount[j + (ncur) * *nClass]== popt1) nodeStatus[ncur]=-1;
      if (tclassCount[j + (ncur +1) * *nClass]== popt2) nodeStatus[ncur+1]=-1;
    }
    
    //on stock les noeuds divises gauche et droite dans treemap
    treemap[0 + (kbuild - 1)  * 2]= ncur +1;
    treemap[1 + (kbuild - 1) * 2]= ncur +2;
    //le noeud kbuild a ete divisé
    nodeStatus[kbuild - 1]=1;
    //on continue pr les autres noeuds
    ncur = ncur+2;
    if (ncur >= *nrnodes) break;   //quitte la boucle 30

  } //fin boucle kbuild
  
  //a partir d ici tous les noeuds ont ete divises, ou sont terminaux (l arbre est construit)
  //ds nbigtree on a le nb noeuds qui existent ds l arbre
  *ndbigtree= *nrnodes;
  for (k = *nrnodes; k>=1; --k){
    if (nodeStatus[k - 1] == 0) *ndbigtree -= 1;
    //on dit que les noeud en attente d etre divises sont teminaux 
    if (nodeStatus[k - 1] ==2) nodeStatus[k - 1] = -1;
  }
  
  //form prediction in terminal nodes
  // nbigtree: nb de noeuds ds l arbre
  for (kn =1; kn <= *ndbigtree; ++kn){
    //si le noeud est terminal
    if (nodeStatus[kn - 1] ==-1){
      //on compte le nb d obs PAR classe (pp) et on assigne au noeud kn la classe ou pp est max
      pp=0;
      for (j=1; j<= *nClass; ++j){
	if (tclassCount[j - 1 + (kn - 1) * *nClass] >pp) {
	  nodeclass[kn - 1]=j;
	  pp = tclassCount[j - 1 + (kn - 1) * *nClass];
	} //end if
	//Break ties at random: s il y a des pp egaux 
	if (tclassCount[j - 1 + (kn - 1) * *nClass] == pp){
	  tieRand = unif_rand();
	  if (tieRand > 0.5){
	    nodeclass[kn - 1]=j;
	    pp=tclassCount[j - 1 + (kn - 1)* *nClass];
	  }//end if
	} //end if tieRand
      } //end for j
    } //end if nodeStatus
  } // for kn

}


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
		   ) {
  
  /*
    commentaires:
    For the best split, msplit is the variable split on. decsplit is the
    dec. in impurity.  If msplit is numerical, nsplit is the case number
    of value of msplit split on, and nsplitnext is the case number of the
    next larger value of msplit.  If msplit is categorical, then nsplit is
    the coding into an integer of the categories going left.
  */

  int last, lcat, mvar, nc,m, i, j, k, l,n;
  double parentNum, parentDen,crit0,critmax, rightNum, rightDen, leftNum, leftDen, u,crit;

  /* compute initial values of numerator and denominator of Gini */
  parentNum = 0.0; // la somme du carré du nombre d'individu/classe
  parentDen = 0.0; // pareil sans le carré
  for (n = 0; n < *nClass; ++n) {
    parentNum += classCount[n] * classCount[n];
    parentDen += classCount[n];
  }
  crit0 = parentNum / parentDen;
  *splitStatus = 0;
  critmax = -1.0e25; 

  
  /* start main loop through variables to find best split. */
  for (j = 0; j < *mtry; ++j) {
    //find the best split on the variables omega
    mvar=omega[j];
    //printf("mvar  %d \n",mvar);

    /* Split on a numerical predictor. */
    rightNum = parentNum;
    rightDen = parentDen;
    leftNum = 0.0;
    leftDen = 0.0;
    zeroDouble(wl, *nClass); //compteur du nombre d'individu par classe à gauche
    for (k = 0; k < *nClass; ++k) wr[k] = classCount[k]; //idem à droite
    for (l = *ndStart; l <= *ndEnd - 1; ++l) { // nd.. permette de connaitre pour le noeud les individus concernés
      nc = a[mvar + (l - 1) * *mDim];
      u = weight[nc - 1]; // récupération de la pondération
      m = class[nc - 1]; // récupération de la classe
      leftNum += u * (2 * wl[m-1] + u);
      rightNum += u * (-2 * wr[m-1] + u);
      leftDen += u;
      rightDen -= u;
      wl[m-1] += u;
      wr[m-1] -= u;
      // on ne s'intéresse uniquement qu'à un changement de facteur
      if (b[mvar + (nc-1) * *mDim] < b[mvar + (a[mvar + l * *mDim]-1) * *mDim]) {
	if (fmin2(rightDen, leftDen) > 1.0e-5) {
	  crit = (leftNum / leftDen) + (rightNum / rightDen);
	  if (crit > critmax) {
	    *nbest = l;
	    critmax = crit;
	    *splitVar = mvar + 1; 
	  }
	  /* Break ties at random: */
	  if (crit == critmax && unif_rand() > 0.5) {
	    *nbest = l;
	    critmax = crit;
	    *splitVar = mvar + 1;
	  }
	}
      }
    }//fin l
    
  } //fin mtry
  if (critmax < -1.0e10 || *splitVar == 0) {
    *splitStatus = -1;
  }
  *decGini = critmax - crit0;
  //printf("bestsplit %d\n", nbest);
}
//#endif /* C_CLASSTREE */

void movedata(int *a, 
	      int *ta,
	      int *mDim,  //nb de variables
	      int *nsample, //nb d obs
	      int *ndStart, 
	      int *ndEnd,
	      int *idmove, 
	      int *ncase, 
	      int *splitVar, 
	      int *nbest, 
	      int *ndendl) { 

/*     a: au debut tableau avec variables en ligne rangees ds ordre croissant et num de l indiv correspondant
c     ta: vecteur de 'transfert' de chaque ligne de a
c     ndstart: indice de depart sur colonnes de a
c     ndend: indice de fin sur colonnes de a
c     idmove: =1 si obs va a gch et 0 si a droite
c     splitVar: renvoye par findbestsplit, meilleure variable pr le split
c     nbest: indice sur tableau a ou a lieu la division du noeud 
      ndendl me semble pas tres utile ici
  */

  int nsp, nc, msh, k, n, ih;
  //initialisation
  zeroInt(idmove, *nsample);
  zeroInt(ta, *nsample); 

  //avant nbest, les obs iront a gauche idmove=1
  for (nsp = *ndStart; nsp <= *nbest; ++nsp) {
    nc= a[(*splitVar-1) + (nsp -1)* *mDim];  //nc est le num d'1 obs
    idmove[nc - 1 ] = 1;    //idmove cree en C
  }
  //apres nbest, les obs iront a droite idmove=0
  for (nsp = *nbest + 1; nsp <= *ndEnd; ++nsp) {
    nc= a[(*splitVar-1) + (nsp -1)* *mDim];  //nc est le num d'1 obs
      idmove[nc - 1] = 0;    //idmove cree en C
  }
  *ndendl=*nbest;
  //printf("%d\n", *ndendl);

  //on change de place les obs qui vont a gauche de celles qui vont a droite
  //on parcourt chaque variable msh
  for (msh = 0; msh <= *mDim-1; ++msh) {
    k = *ndStart -1;
      // on parcourt les lignes de a
      for (n = *ndStart; n <= *ndEnd; ++n){
	ih = a[msh + (n-1)* *mDim];  // ih num d un indiv
	// les elements ayant idmove=1 sont d abord stocke ds ta
	if ( idmove[ih - 1] ==1){
	  k += 1;  //k commence a ndStart
	  //ta recoit le num des indiv a gauche
	  ta[k]=a[msh + (n-1)* *mDim];
	}
      } // end for n
      for (n = *ndStart; n <= *ndEnd; ++n){
	 ih = a[msh + (n-1)* *mDim];
	// on stocke ensuite ds ta les elements ayant idmove=0
	if ( idmove[ih - 1] ==0){
	  k += 1;
	  ta[k]=a[msh + (n-1)* *mDim];
	}
      } // end for n
      // les elements ayant idmove=1 sont rerange a gche de a puis les elements idmove=0 a droite, et ce pour CHAQUE variable
      for ( k = *ndStart; k <= *ndEnd; ++k){
	a[msh + (k-1)* *mDim] = ta[k];
      }
  }  // end for msh


  //compute case nos. for right and left nodes
  // cas de splitVar continue:
  // on met ncase a jour pr la variable splitVar
  for (n = *ndStart; n <= *ndEnd; ++n){
    ncase[n-1]= a[(*splitVar-1) + (n-1)* *mDim];
  } 

 
}  // fin movedata


