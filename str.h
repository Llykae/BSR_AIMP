/*--------------------------------------------------------------------------
   This file contains all global variables for the Pseudo-potential code.
--------------------------------------------------------------------------*/
#ifndef STR_H
#define STR_H


#define BSLMAX 10                               //twice largest angular momentum in basis set + 4 -- d ==> BSLMAX = 8
#define BSLMAX2 2*BSLMAX
#define BSLMAX3 3*BSLMAX
#define BSLMAX6 6*BSLMAX
#define NBSFMAX 256                             //max number of primitive basis functions per atoms
#define NBSFMAX1 NBSFMAX+1
#define NBCCMAX 50                              //max number of radial basis functions for contraction
#define NBCCMAX1 NBCCMAX+1

struct gauss_basis_set{
  double ne;                                    //number of electron
  int nfs, nfp, nfd;				//number of function for orbital 1s, 2p, ...
  int ngs, ngp, ngd;				//number of gaussian for 1s, 2p, ...
};

struct orbital{
  int l;                                        //angular momentum
  double arg[NBCCMAX];                          //exponant
  double cgauss[NBCCMAX];                       //gaussian coeff
  double nrm[NBCCMAX];                          //gaussian norm
  double cdensity[NBCCMAX][NBCCMAX];            //coefficient overlap for density matrix
  double cexchange[NBCCMAX][NBCCMAX];           //coefficient overlap for exchange matrix
};

struct bsp
{
  double BspValue[48];
};

struct pgl
{
  //double pvs[8][24][48];
  double pvs[12][48];
  double pvp[12][48];
  double pvd[12][48];
}; 

struct phase_shift
{
 double ps[1000];
};

struct atom
{
  double occ[6];
};

struct basis_extraction
{
  int size_be;
  int ns,np,nd;
  double arg[NBCCMAX];
  double args[NBCCMAX],argp[NBCCMAX],argd[NBCCMAX];
  double nrm[NBCCMAX];
  double vs[NBCCMAX][NBCCMAX];
  double vp[NBCCMAX][NBCCMAX];
  double vd[NBCCMAX][NBCCMAX];
};

struct vshell					//Shell structure for AIM Pseudo Potential
{
   int l;
   double a;
   double norm;
};

struct vxaimp
{
   double zi;					//the ionic charge
   double qe;					//the core electronic charge
   double zc;					//net core charge : zc = zi - qe
   double epauli;
//---
   int nlb[4];
   int lbmax;
   struct vshell *shlb[4];
   double **rho[4];				//at most lmax = 3 for the core
//---
   int nlx[6];
   int lxmax;
   struct vshell *shlx[6];
   double **vvx[6];				//at most lmax = 5 for the core
};

struct potcor
{
    int size_be;
    int lmax;
    double arg[20];
    double vab[20][20];
};


struct param_pola
{
    double aab[5];        //alphad, alphaq, beta
    char *atom;

    int ndipol[3];
    int *dipow[3];
    double *dipc[3];
    double *dipg[3];

    int nquad[4];
    int *quadpow[4];
    double *quadc[4];
    double *quadg[4];

    int noctu;
    int *octupow;
    double *octuc;
    double *octug;

};

#endif

