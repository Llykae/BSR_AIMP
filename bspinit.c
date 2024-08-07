/*--------------------------------------------------------------
   B-spline and atomic data initialization:

		B. Gervais --- 2012
		B. Gervais + E. Hochard --- 2021
--------------------------------------------------------------*/
#define NKNOT_MAX 10000				//no more than 1000 points expected

#include "bsputil.h"

void init_data(double *z, int *ne, double *q, int icontinuum);
void init_basis(int icontinuum, int kb, int nspline, double *tn, int *nbasis, struct bsp_basis *bss);
void mk_knot(double zat);
void mk_knot_cont(double zat);
void read_bs_size(char in_knts[256], int *kb, int *nbp, int *nb);
void read_knot(char in_knts[256], int nb, int kb, int nbp, double *tn, double *rn);



/*--------------------------------------------------------------
   Atomic data initialization
--------------------------------------------------------------*/
void init_data(double *z, int *ne, double *q, int icontinuum)
{ 
   *q = *z - (double)(*ne); 

   if(icontinuum == 0) mk_knot(*z);
   else mk_knot_cont(*z);
}



/*--------------------------------------------------------------
   Basis set initialisation. There ae 2 main options defined by
   the boundary conditions:

   1-icontinuum = 0:	f(0)		= 0.0
			f(xmax)		= 0.0
   In such a case, the basis functions are identical to the
   B-splines, but:
   -the first B-splines b1(x), for which b1(x=0) = 1 is removed
   from the basis set.
   -the last k splines, where k is the degree of the slpines,
   are also removed, so that only the last Bsplines contributes
   on the boundary, with a value 0.0, due to (x-xmax)^k behaviour
   on the limit xmax.

   2-icontinuum = 1:	f(0)		= 0.0
			(f'/f)(xmax)	!= 0.0
   In such a case, the basis functions are identical to the
   B-splines, but:
   -the first B-splines b1(x), for which b1(x=0) = 1 is removed
   from the basis set.
   the last k splines are  kept and the basis functions verify:
	fN(xmax) = BN-1(xmax) = 1.0
	dfN/dx(xmax) = dBN-1/dx(xmax) = +k

	fN-1(xmax) = BN-1(xmax) = 0.0
	dfN-1/dx(xmax) = dBN-1/dx(xmax) = -k

	Bm(xmax) = Bm'(xmax) = 0.0	for m < N-1
--------------------------------------------------------------*/
void init_basis(int icontinuum, int kb, int nspline, double *tn, int *nbasis, struct bsp_basis *bss)
{
   int m,n;
   //double bmin,bmax,dbmin,dbmax;

   switch(icontinuum)
      {
      case 0:						//generates nbsplines-1 functions fn with fn(xmax) = fn(0) = 0.0
         printf("icontinuum = %d \n",icontinuum);
	 for(n=2;n<nspline;n++)
            {
            m = n - 1;
            bss[m].ib = n;				//basis n-1 is bsplines n
            bss[m].tmin = tn[n];
            bss[m].tmax = tn[n+kb+1];
            /*bmin = bspline(tn,n,kb,bss[m].tmin);
            bmax = bspline(tn,n,kb,bss[m].tmax);
            dbmin = dbspline(tn,n,kb,bss[m].tmin);
            dbmax = dbspline(tn,n,kb,bss[m].tmax);*/
//printf("bss: ib=%d \t tmin=%+le \t tmax=%+le \t bmin=%+le \t bmax=%+le \t dbmin=%+le \t dbmax=%+le \n",bss[n-1].ib,bss[n-1].tmin,bss[n-1].tmax,bmin,bmax,dbmin,dbmax);
            }
         *nbasis = nspline - 2;	
         break;

      case 1:
         printf("icontinuum = %d \n",icontinuum);
	 for(n=2;n<nspline-1;n++)
            {
            m = n - 1;
            bss[m].ib = n;				//basis n-1 is bsplines n
            bss[m].tmin = tn[n];
            bss[m].tmax = tn[n+kb+1];
            /*bmin = bspline(tn,n,kb,bss[m].tmin);
            bmax = bspline(tn,n,kb,bss[m].tmax);
            dbmin = dbspline(tn,n,kb,bss[m].tmin);
            dbmax = dbspline(tn,n,kb,bss[m].tmax);*/
//printf("bss: ib=%d \t tmin=%+le \t tmax=%+le \t bmin=%+le \t bmax=%+le \t dbmin=%+le \t dbmax=%+le \n",bss[n-1].ib,bss[n-1].tmin,bss[n-1].tmax,bmin,bmax,dbmin,dbmax);
            }
         *nbasis = nspline - 3;	
         break;

      default:
         printf("no such possibility: icontinuum = %d \n",icontinuum);
         exit(1);
         break;
      }
}



/*--------------------------------------------------------------
   B-spline generation. The output is saved in the file
   "bspline_knot.dat" to be read further.

   Original parameters from Zatsarinny and Froese-Fischer
   (J. Phys. B 33 p.313-341) good for 2 S-states of H:
	mb = 3
	kb = 8
	hmax = 1.0
	tmax = 25.0*zat
	dtmax = hmax = 1.0
--------------------------------------------------------------*/
void mk_knot(double zat)
{
   int i,imax;
   int mb;
   int kb;
   int nbp;
   double h,hmax;
   double dt,dtmax,tmax;
   double tn[NKNOT_MAX+1];
   double rn[NKNOT_MAX+1];
   FILE *fout;

   mb = 3;
   kb = 8;
   h = powint(2.0,-mb);
   hmax = 1.0;
   tmax = 100.0*zat;
   dtmax = 2.0;

   for(i=1;i<=kb;i++) tn[i] = 0.0;
   while(tn[i] < hmax)
      {
      tn[i+1] = tn[i] + h;
      i++;
      }
   dt = 0.0;
   while(dt < dtmax)
      {
      tn[i+1] = tn[i]*(1.0 + h);
      dt = tn[i+1] - tn[i];
      i++;
      }
   while(tn[i] < tmax)
      {
      tn[i+1] = tn[i] + dtmax;
      i++;
      }
   imax = i;
   nbp = imax;

   fout = fopen("bspline_knot.dat","w");
   fprintf(fout,"%d \t\t\t\t B-spline order \n",kb);
   fprintf(fout,"%d \t\t\t\t number of knots \n",nbp);
   for(i=1;i<=imax;i++)
      {
      rn[i] = tn[i]/zat;
      fprintf(fout,"%1.12le \t %1.12le \n",tn[i],rn[i]);
      }
   fclose(fout);

   if(nbp >= NKNOT_MAX)
      {
      printf("number of knots exceeds %d --- exit \n",NKNOT_MAX);
      exit(1);
      }
}



/*--------------------------------------------------------------
   Basis spline generation for continuum function generation in
   elastic scattering computations. The output is saved in the
   file "bspline_knot.dat" to be read further.

   Original parameters from Zatsarinny and Froese-Fischer
   (J. Phys. B 33 p.313-341) good for 2 S-states of H in the inner
   part + modification for continuum.
	mb = 3
	kb = 8
	hmax = 1.0
	tmax = pi*zat*(3 + E(zat^1/3))
	dtmax = hmax = 1.0
   NB: t rescaled to give r = t/zat
--------------------------------------------------------------*/
void mk_knot_cont(double zat)
{
   int i,j,imax;
   int mb;
   int kb;
   int nbp;
   double h,hmax;
   double dt,dtmax,tmax;
   double tn[NKNOT_MAX+1];
   double rn[NKNOT_MAX+1];
   FILE *fout;

   mb = 3;
   kb = 8;

   h = powint(2.0,-mb);
   hmax = 1.0;
   tmax = 4.0*zat;
   //tmax = 10.0*acos(-1.0)*zat*(3.0 + ceil(pow(zat,1.0/3.0)));
   //tmax = 1.1*acos(-1.0)*zat*(3.0 + ceil(pow(zat,1.0/3.0)));
   dtmax=1.0;

   for(i=1;i<=kb;i++) tn[i] = 0.0;
   while(tn[i] < hmax)
      {
      tn[i+1] = tn[i] + h;
      i++;
      }
   dt = 0.0;
   while(dt < dtmax)
      {
      tn[i+1] = tn[i]*(1.0 + h);
      dt = tn[i+1] - tn[i];
      /*printf("%lf %lf\n",tn[i+1],dt);
      getchar();*/
      i++;
      }
   while(tn[i] < tmax)
      {
      tn[i+1] = tn[i] + dtmax;
      /*printf("%lf %lf\n",tn[i],tmax);
      getchar();*/
      i++;
      }
//---
//   tn[i+3] = tn[i+2] = tn[i+1] = tn[i];			//add kb more knots to finish
//   i+=3;
   for(j=1;j<=kb+1;j++)						//add kb+1 more knots to finish
      {
      tn[i+1] = tn[i];
      i++;
      }
   imax = i;
   nbp = imax;

   fout = fopen("bspline_knot.dat","w");
   fprintf(fout,"%d \t\t\t\t B-spline order \n",kb);
   fprintf(fout,"%d \t\t\t\t number of knots \n",nbp);
   for(i=1;i<=imax;i++)
      {
      rn[i] = tn[i]/zat;
      fprintf(fout,"%1.12le \t %1.12le \n",tn[i],rn[i]);
      }
   fclose(fout);

   if(nbp >= NKNOT_MAX)      {
      printf("number of knots exceeds %d --- exit \n",NKNOT_MAX);
      exit(1);
      }
}




/*--------------------------------------------------------------
   Basis splines generation. The number of knots and the degrees
   of the splines is read in the input file in_knts.
--------------------------------------------------------------*/
void read_bs_size(char in_knts[256], int *kb, int *nbp, int *nb)
{
   FILE *fin;
   char strtmp[1024];

   fin = fopen(in_knts,"r");
   fscanf(fin,"%d\t",kb);
   fgets(strtmp,256,fin);
   fscanf(fin,"%d\t",nbp);
   fgets(strtmp,256,fin);

   *nb = *nbp - *kb;
   printf("basis set parameters: nb=%d \t kb=%d \t nbp=%d \n",*nb,*kb,*nbp);

   fclose(fin);
}



/*--------------------------------------------------------------
   Basis splines generation. The series of knots is read in the
   input file in_knts and stored in arrays tn and rn.
--------------------------------------------------------------*/
void read_knot(char in_knts[256], int nb, int kb, int nbp, double *tn, double *rn)
{
   FILE *fin;
   char strtmp[1024];
   int i;

   fin = fopen(in_knts,"r");
   fgets(strtmp,256,fin);
   fgets(strtmp,256,fin);

   for(i=1;i<=nbp;i++)
      {
      fscanf(fin,"%le \t %le \n",&tn[i],&rn[i]);
      //printf("i=%d \t ti=%1.12e \t ri=%1.12e \n",i,tn[i],rn[i]);
      }

   fclose(fin);
}
