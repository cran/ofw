#
# Copyright (C) 2007 Kim-Anh Lê Cao, Patrick Chabrier, INRA,
# French National Institute for Agricultural Research.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#this function projects all elements in P on a simplex

simplexe=function(P)
	{
	#1. Project on hyperplane Hf if P does not belong to Hf
	somme=sum(P)
	lambda=(somme-1)/length(P)
	P=P-lambda  #projection on Hf
	
	#2. Project on the simple the elements of P <0
	Sf= 1		#if not on the simplex
	while(Sf!=0)
	{
		Sf=0
		J=sum(P<=0)          #nb of P elements <=0
		somme=sum(P[P>0])    
		P[P<=0]=0            #negative elements are set to zero
		P[P>0]=P[P>0] + (1-somme)/(length(P) - J)   #projection on the simplex
		Sf=sum(P<0)    #are there any other elements not on the simplex?
	}
	return(P)
}
