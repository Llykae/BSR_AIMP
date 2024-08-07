/*---------------------------------------------------------------------
	Main program for potential extraction from atomic orbitals.
	The code extract:
		The density matrix in orbital basis
		A Pauli repulsion potential in orbital basis
		An exchange potential in additional extraction basis
-----------------------------------------------------------------------
	The code is made for a single center potential on one single
	atom. The extraction of exchange is thus made in primitive
	spherical hamonics with radial GTO. The Coulomb integrals are
	computed accordingly.
-----------------------------------------------------------------------
	Atomic units are used throughout. Implicitely we thus have:
		electron mass:		me = 1
		planck constant:	h/2pi = 1
		light velocity:		c = 137
		nucleon mass:		mn = 1826
		Boltzmann constant:	kb = 3.166089e-06
-----------------------------------------------------------------------
	B. GERVAIS and E. HOCHARD  12/2021
---------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

#include "genaimp.h"						//global variables and structures for extraction code 
#include "globmath.c"						//global mathematical variables for LCAO code
#include "nrutil.c"						//dynamical memory allocation functions

#include "pfmath.c"						//some mathematical functions
#include "diag.c"						//diagonalisation functions

#include "genaimpinit.c"					//initialisation file
//---------------------

#include "extraction.c"

#define EPAULI 300.0



int main()
{
//-----------physical variables-----------------
   int natm;							//number of atoms in the cluster

   int nbct;
   int nbctcs;							//total number of contracted spherical basis functions
   int nbct0;							//total number of primary cartesian basis functions
   int lbmax;
   int lxmax;
   double zat,qat;						//atomic charges (nucleus, effective)

//-----------data variables---------------------
   struct atom_dat *at_dt;					// the structure of atom data
   struct shell_basis *shbss[LBASMAX];
   struct shell_basis *shx[LEXTMAX];				// Extraction shell
   int nlshx[LEXTMAX];
   int nlshb[LBASMAX];						//at most l = LBASMAX - 1 - l in [0,LBASMAX-1]
   double **rho0[LBASMAX];					//1-body density matrix in atomic basis for each l value

//-----------output variables-------------------	
   FILE *fout;

//-----------other variables----------------
   int i,j,k;
   int lambda;
   double sum; 
   double epauli = EPAULI;


   initmath();									//mathematical constant initialization
   natm = 1;

   printf("----------------------------\n");
   printf("\t Number of atoms : natm = %d \n",natm);

//------------ALLOCATE MEMORY-------------------------
   at_dt=malloc((natm+1)*sizeof(struct atom_dat));  				//allocate memory for atom data structure at_dt
   initatm(natm,&at_dt);							//get atom data
   readorbat(nlshb,shbss,&lbmax,&zat,&qat,rho0,at_dt[1].bsnm);			//get orbitals and associated basis shell

//-------------------
   fout = fopen("output/vaimp.out","w");					//write density matrix for Pauli repulsion evaluation
   fprintf(fout,"# AIMP for closed-shell spherically symmetric sytem\n");
   fprintf(fout,"# Nuclear charge - effective charge\n");
   fprintf(fout,"ATOM\n");
   fprintf(fout,"%1.2lf\t%1.2lf\n",zat,qat);
   fprintf(fout,"# Density matrix from input occupied orbitals\n");
   fprintf(fout,"DENSITY\n");
   fprintf(fout,"%d\n",lbmax);
   for(lambda=0;lambda<=lbmax;lambda++)
      {
      fprintf(fout,"# l\tnl\n");
      fprintf(fout,"%d\t%d\n",lambda,nlshb[lambda]);
      for(i=1;i<=nlshb[lambda];i++) fprintf(fout,"%+1.4le  \t",shbss[lambda][i].a);
      fprintf(fout,"\n");
      for(i=1;i<=nlshb[lambda];i++)
         {
         for(j=1;j<=nlshb[lambda];j++) fprintf(fout,"%+le\t",rho0[lambda][i][j]);
         fprintf(fout,"\n");
         }
      }
   fprintf(fout,"#-----------------------------------------\n");
   fprintf(fout,"# ---------------- vpauli - Ep=%+1.2f a.u.:\n",epauli);
   fprintf(fout,"EPAULI\n");
   fprintf(fout,"%+le\n",epauli);
   fprintf(fout,"#-----------------------------------------\n");
   fclose(fout);
//exit(1);

//---------------------------------------------------
   initxshell(&lxmax,nlshx,shx,"basis/extraction.dat");
   extraction(lxmax,nlshx,shx,lbmax,nlshb,shbss,rho0);

//---------------------------------------------------
   free(at_dt);
   exit(1);
}
