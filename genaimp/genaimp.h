/*--------------------------------------------------------------------
   This file contains all global variables for the LCAO code.
-------------------------------------------------------------------*/
//#define IOPENMP 1				//switch for openomp parallel options

#ifndef GENAIMP_H
#define GENAIMP_H 

//#define MAXPSP 7				//largest number of elements (including their various representations)
//#define MAXELEC 21				//largest number of electrons + 1

#define BSLMAX 8				//twice largest angular momentum in basis set + 4 -- d ==> BSLMAX = 8
#define BSLMAX2 2*BSLMAX
#define BSLMAX3 3*BSLMAX
#define BSLMAX6 6*BSLMAX
#define NBSFMAX 100				//max number of primitive basis functions per atoms
#define NBSFMAX1 NBSFMAX+1
#define NBCCMAX 16				//max number of radial basis functions for contraction
#define NBCCMAX1 NBCCMAX+1
//#define NPROJMAX 180		/*max number of projectors*/
//#define NPROJMAX1 NPROJMAX+1
//#define NPROJMAX1 NPROJMAX+1
//#define NANG_MAX 4				//max number of angular points
//#define NRAD_MAX 256				//max number of radial points
//#define NDS_MAX 10		/*max size of degenerate space*/
//#define NDS_MAX1 NDS_MAX+1

#include "constant.h"				//physical constant and conversion factors

struct atom_dat
{
   double  r[4];				// position vector of the atom
   double  p[4];				// momentum vector of the atom
//   double  f[4];				//the force experienced by the atom
//   double  f0[4];				//previous force (for Beeman)
//   double  fm[4];				//previous force (for Beeman)
   double  d[4];				//
   double  m;					// mass of the atom
   double  qat;					//atomic charge for electric field generation
//   double  alpha;				//dipole polarizability
//   double  beta;				//dynamical correction to dipole polarizability
//   double  qpol;				//quadrupole polarizability
   double  ge;					//electronic cutoff
   double  gi;					//ionic cutoff
   int     zat;					//atomic number index
//   int     itype;				//switch for cluster / matrix
//   int     ibasis;				//switch for basis / no basis 
//   int     idpol;				//switch for dipole active / inactive
//   int     iqpol;				//switch for quadrupole active / inactive
//   int     idpoldpol;				//switch for dipole-dipole  active / inactive
   int     imol;				//Number of the molecule the element belongs to 
//   int     ips;					//index of potential parameters structure
   char    element[3];				//Name (H,He,etc) of the element
   char    bsnm[64];				//atom basis set file name
   char    pspnm[64];				//atom potential file name
};


/*---
//-------------------------------------------------------------------
struct pxyz_basis_set				// structure for cartesian primitive (uncontracted) basis set
{
   int iat;					//unique atom center
   double r[4];					//position vector of the basis function center
   double nrm;					//norm of the basis function
   double a;					//exponent
   int p0;					//l power
   int p1;					//x power
   int p2;					//y power
   int p3;					//z power
   double intg;					//integral of the function (for fitting only)
};


//-------------------------------------------------------------------


//-------------------------------------------------------------------
struct clm_basis_set				// structure for real spherical harmonics contracted basis set
{
   int iat;					//unique atom center
   int ir;					//index of radial atomic basis function
   int l;					//l angular momentum = spherical harmonic l index
   int m;					//m real spherical harmonic m index
   int nmxyz;					//number of cartesian expansion coefficients (max = BSLMAX2)
   int indx[BSLMAX2];				//index of cartesian expansion coefficients (max = BSLMAX2)
   double cmxyz[BSLMAX2];			//coefficients of cartesian expansion coefficients (max = BSLMAX2) 
   double nrm;
};


//-------------------------------------------------------------------
struct basis_set
{
   int iat;					//unique atom center
   double r[4];                 //position vector of the basis function center
   double nrm;                  //norm of the basis function
   double a;                    //exponent
   int p0;                      //l power
   int p1;                      //x power
   int p2;                      //y power
   int p3;                      //z power
};
---*/


//-------------------------------------------------------------------
struct shell_basis				//shell structure
{
   int iat;					//atomic index
   int l;					//angular momentum
   double a;					//exponent - for cartesian functions
   int ir;					//radial function index for spherically contracted functions
   int n;					//number of functions in the shell : n = (l+1)*(l+2)/2 for cartesian GTO
   int *i;					//index of basis function in the shell
};

#endif
                
