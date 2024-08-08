/*--------------------------------------------------------------
   A B-spline programm for Atomic calculations

			B. Gervais- E. Hochard --- 2021
--------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"

#include "globmath.c"
#include "gauss.c"
#include "bspmath.c"
#include "diag.c"
#include "bessel.c"

#define ICONTINUUM 1			                    //set continuum switch to ensure outer boundary condition
#define NOBLE_ATOM 0                          //0 = He / 1 = Ne / 2 = Ar / 3 = Kr / 4 = Xe
#define LANG 0						                    //ANGULAR MOMENTUM
#define LMAX 2
#define Nk 500                                //NUMBER OF K POINT FOR TCS EVALUATION (500 = 13.605eV APPROXIMATIVELY)
#define DKH 0                                 //SIGN = 0, Without DKH. SIGN = 1, With DKH
#define TDSWC 0                               //0 = TCS  -   1 = DCS

#define SIZE 2                                //NUMBER OF K POINT (IN TERMS OF ENERGY)
#define DCS_LMAX 3                            //NUMBER OF ANGULAR MOMENTUM IN THE INPUT FILE - CAN'T BE CHANGED FOR NOW
#define EXTRACTION_DCS 1                      //0 FOR YES - 1 FOR NO
#define EXTRACTION_DCS_NAME "Williams"        //OUTPUT LOCATION FOR DCS EXTRACTION - ONLY NAME - MAKE SURE TO CREATE DIRECTORY BEFORE
#define EXTRACTION_DCS_FILE "Williams.data"                //FILE WITH PHASE-SHIFT TO USE FOR DCS EVALUATION

#define POLA  0                               //IN CASE YOU WANT DATA PLOT FOR POLARISATION POTENTIAL CHANGE FOR 1


#include "str.h"
#include "readbs.c"
#include "readbe.c"

#include "gto.c"
#include "bspinit.c"
#include "bsputil.c"
#include "hmlt_function.c"
#include "Polarisation.c"
#include "bsphmlt.c"
#include "GL_function.c"

#include "locpotential.c"
#include "correlation_function.c"

#include "Renorm.c"
#include "bspvex.c"
#include "bspvexact.c"

#include "Legendre_Polynomials.c"
#include "bsprelat.c"
#include "DKH_BSP.c"

#include "calogero.c"
#include "scattering_function.c"
#include "main_function.c"



int main()
{

  int icontinuum=ICONTINUUM;

  int ne;						                      //number of electron
  double zat, qat;						            //atomic number (nuclear charge) + atomic charge

  int l,swc,i;
  int nbasis, nspline;						        //number of basis functions + number of splines - can be different from nbasis (nb = ns-1)
 
  int kb, nbp;						                //degree of B-splines + number of knots
  double *tn, *rn;						            //the array of knots + the array of corresponding radius

  double E[Nk];
  
  double **ss,**tt,**vvx,**vvr, **vvc,**hh1;
  double **coulomb, **ionique, **pola;
  double *eorb,**corb;

  char *fname;
  char *data;

  struct bsp_basis *bss;
  struct gauss_basis_set *gbs;
  struct orbital *orb;
  struct basis_extraction *be;
  struct bsp *bspline;
  struct atom *at;
  struct phase_shift *d0;

  struct vxaimp pspa;
  struct potcor *vcc;
  struct param_pola *ppola;

  double **dkh_hdk2, **vcrelat;

 

  ppola=malloc(sizeof(struct param_pola));
  switch(NOBLE_ATOM)
  {
    case 0 :
      ppola->atom = "He";
    break;

    case 1 :
      ppola->atom = "Ne";
    break;

    case 2 :
      ppola->atom = "Ar";
    break;

    case 3 :
      ppola->atom = "Kr";
    break;
      
    case 4 :
      ppola->atom = "Xe";
    break;      
  }
  
  printf("\n                              START READING BASIS AND POTENTIAL\n\n");
  readvaipp(&pspa,ppola);
  printf("\n\n\n                                 START PHASE SHIFT ANALYSIS\n\n");

  data = (char *) malloc(200*sizeof(char));
  sprintf(data,"./basis_set/%s_basis_set.dat",ppola->atom);

  gbs = malloc(sizeof(struct gauss_basis_set));  
  init(gbs,data,ppola->atom);

  zat=gbs[0].ne;
  ne=zat;
  init_gauss();
  init_data(&zat,&ne,&qat,icontinuum);
  printf("atomic parameters: zat=%lf \t ne=%d \t qat=%lf \n",zat,ne,qat);



  read_bs_size("bspline_knot.dat",&kb,&nbp,&nspline);
  tn = dvector(0,nbp);
  rn = dvector(0,nbp);
  read_knot("bspline_knot.dat",nspline,kb,nbp,tn,rn);

  
  bss = malloc((nspline+1)*sizeof(struct bsp_basis));
  init_basis(icontinuum,kb,nspline,tn,&nbasis,bss);



  orb = malloc((gbs[0].nfs+gbs[0].nfp+gbs[0].nfd+1)*sizeof(struct orbital));
  bspline = malloc((nspline+1) * sizeof(struct bsp));
  be=malloc(sizeof(struct basis_extraction));
  at=malloc(sizeof(struct atom));
  d0=malloc((LMAX+1)*sizeof(struct phase_shift));
  vcc=malloc((LMAX+1)*sizeof(struct potcor));

  swc=0;
  if((gbs[0].ne>4) && (gbs[0].ne<=18)) {swc=1;}
  else if(gbs[0].ne>18)                {swc=2;}

  
  ss = dmatrix(1,nbasis,1,nbasis);			///allocate overlap matrix
  tt = dmatrix(1,nbasis,1,nbasis);			//allocate kinetic energy matrix
  pola = dmatrix(1,nbasis,1,nbasis);
  coulomb = dmatrix(1,nbasis,1,nbasis);
  ionique = dmatrix(1,nbasis,1,nbasis);
  vvx = dmatrix(1,nbasis,1,nbasis);                     //exchange potential
  vvr = dmatrix(1,nbasis,1,nbasis);                     //pauli potential
  vvc = dmatrix(1,nbasis,1,nbasis);
  hh1 = dmatrix(1,nbasis,1,nbasis);			//allocate 1-electron hmiltonian matrix

  eorb = dvector(1,nbasis);
  corb = dmatrix(1,nbasis,1,nbasis);			//allocate orbital coefficient matrix

  dkh_hdk2 = dmatrix(1,nbasis,1,nbasis);
  vcrelat = dmatrix(1,nbasis,1,nbasis);
  
  read_input_pola(ppola);

  for(l=LANG;l<=LMAX;l++){ 
    build_potential(zat,nbasis,nspline,kb,l,tn,ss,tt,pola,ionique,coulomb,vvx,vvr,vvc,bss,bspline,be,&pspa,vcc,ppola,dkh_hdk2,vcrelat);
    build_eigen(l,nbasis,eorb,corb,ss,tt,pola,ionique,coulomb,vvx,vvr,vvc,hh1,dkh_hdk2,vcrelat);
    build_PS(l,nbasis,E,eorb,corb,d0,bss,&pspa,ppola);}
  if(TDSWC==0) build_TCS(E,d0);
  else build_DCS(ppola,d0);

  if(EXTRACTION_DCS == 0) build_DCS_bis(ppola);
  return 0;
}
