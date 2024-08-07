/*----------------------------------------------------------------------------
   Phase shift analysis for elastic diffusion of electron in a spherical
   potential. Atomic units are used everywhere. rescaling of the potential
   is done to account for condensed phase, and the potential is set to 0
   at distances larger then rc.
-----------------------------------------------------------------------------*/

//#include "globmath.c"
//#include "bessel.c"						//calculation of spherical Bessel function


//#define M 60							//number of velocities
//#define M1 M+1
//#define LMAX 6							//number of phase shifts
//#define LMAX1 LMAX+1

//#define IOPT 7							//potential option

#define CALOGERO_RMAX 200.0					//outermost radius for phase shift propagation

#define SAFE 0.9						//parameter for Runge-Kutta integrator
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define PGROW  -0.2      


double dcalogero(int lang, double v, double d0, double r0, struct vxaimp *pspa, struct param_pola *ppola);
void rkqs6(double *y, double *x, double htry, double eps, double *yscal, double *hnext, double v, int l, struct vxaimp *pspa, struct param_pola *ppola);
void rkck6(double *y, double h, double *yout, double *yerr, double v, int l, double x, struct vxaimp *pspa, struct param_pola *ppola);

double derivs(double v, int l, double r, double y, struct vxaimp *pspa, struct param_pola *ppola);
double pot(double r, int lang, struct vxaimp *pspa, struct param_pola *ppola);

double dcalogero(int lang, double v, double d0, double r0, struct vxaimp *pspa, struct param_pola *ppola)
{
   double r,dr,rc;
   double d,dscal;
   double eps=1.0e-10;

   r = r0;								//starting radius
   rc = CALOGERO_RMAX;							//outermost radius
   dr = 0.0001;								//initial integration step for adaptative Runge-Kutta
   
   d = d0;								//initial phase shift == 0
//printf("l=%d \t v=%e \t r=%e \t d=%e \t u=%+le \n",lang,v,r,d,pot(r,lang,pspa));
   while(r<rc)
      {
      dscal = d + dr;
      rkqs6(&d,&r,dr,eps,&dscal,&dr,v,lang,pspa,ppola);
//      printf("l=%d \t v=%e \t r=%e \t d=%e \n",l,v,r,d);
      }
//printf("l=%d \t v=%e \t r=%e \t d=%e \t u=%+le \n",lang,v,r,d,pot(r,lang,pspa));
//getchar();

   return(d);
}


/*--------------------------------------------------------------------------------
   rkqs6 defines the current step and gives an estimation of the next step to be
   done to integrate a sdifferential equations by a fifth order Runge-Kutta method
   (see numerical recipies in C). The Runge-Kutta step is done in the routine
   rkck6() which also provide an estimate of the error.
--------------------------------------------------------------------------------*/
void rkqs6(double *y, double *x, double htry, double eps, double *yscal, double *hnext, double v, int l, struct vxaimp *pspa, struct param_pola *ppola)
{
   double errmax,h,htemp,xnew,yerr,ytemp;
  
   h=htry;
   for(;;)
      {
      rkck6(y,h,&ytemp,&yerr,v,l,*x,pspa,ppola);
      errmax = 0.0;
      errmax = MAX(errmax,fabs(yerr/(*yscal)));
      errmax /= eps;
      if(errmax > 1.0)
         {
         htemp = SAFE*h*pow(errmax,PSHRNK);
         h=(h>=0.0 ? MAX(htemp,0.1*h):MIN(htemp,0.1*h));
         xnew=(*x)+h;
         if (xnew == *x) printf("step too small");
         continue;
         }
      else 
         {
         if (errmax>ERRCON) *hnext=SAFE*h*pow(errmax,PGROW);
         else *hnext = 5.0*h;
         *x += h;
         *y = ytemp;
         break;
         }
      }
}
   
void rkck6(double *y, double h, double *yout, double *yerr, double v, int l, double x, struct vxaimp *pspa, struct param_pola *ppola)
{
   double derivs();
   static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,
                 b21=0.2,b31=3.0/40.0,b32=9.0/40.0,
                 b41=0.3,b42=-0.9,b43=1.2,
                 b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
                 b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
                 b64=44275.0/110592.0,b65=253.0/4096.0,
                 c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
                 dc5=-277.0/14336.0;
   double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   double ak1,ak2,ak3,ak4,ak5,ak6;
   double xtemp,ytemp;

   ak1 = derivs(v,l,x,*y,pspa,ppola);							//variation
   ytemp = *y + b21*h*ak1;

   xtemp = x + a2*h;
   ak2 = derivs(v,l,xtemp,ytemp,pspa,ppola);				//variation
   ytemp = *y + h*(b31*ak1 + b32*ak2);

   xtemp = x + a3*h;
   ak3 = derivs(v,l,xtemp,ytemp,pspa,ppola);				//variation
   ytemp = *y + h*(b41*ak1 + b42*ak2 + b43*ak3);

   xtemp = x + a4*h;
   ak4 = derivs(v,l,xtemp,ytemp,pspa,ppola);				//variation
   ytemp = *y + h*(b51*ak1 + b52*ak2 + b53*ak3 + b54*ak4);

   xtemp = x + a5*h;
   ak5 = derivs(v,l,xtemp,ytemp,pspa,ppola);				//variation
   ytemp = *y + h*(b61*ak1 + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5);

   xtemp = x + a6*h;
   ak6 = derivs(v,l,xtemp,ytemp,pspa,ppola);				//variation
   *yout = *y + h*(c1*ak1 + c3*ak3 + c4*ak4 + c6*ak6);

   *yerr = h*(dc1*ak1 + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);
}


double derivs000(double v, int l, double x, double y)
{
   double dydx;

   dydx = exp(-x);

   return(dydx);
}


double derivs(double v, int l, double r, double y, struct vxaimp *pspa, struct param_pola *ppola)
{
   double jl,nl,djl,dnl;
   double t,u,dydx,z;
   int m;

   z = v*r;
   u = 2.0*pot(r,l,pspa,ppola);

   if(z<1.0e-6)								//t = r*jl*cos(y)
      {
      t = r*cos(y);
      for(m=1;m<=l;m++) t *= z/(2*m + 1);
      }
   else
      {
      bslsph(l,z,&jl,&nl,&djl,&dnl);
//printf("r=%e \t rc=%e \t v=%e \t l=%d \n",r,rc,v,l);
//printf("jl=%e \t nl=%e \t djl=%e \t dnl=%e \n",jl,nl,djl,dnl);
      t = r*(jl*cos(y) - nl*sin(y));
      }
   dydx = -v*u*t*t;

   return(dydx);
}



//=====================================================

/*----------------------------------------------------------------------
   calculation of the interaction potential between the projectile
   electron and the target atoms.
	input:	r	electron-nucleus distance
		lang	angular momentum
	output:	v	the potential value
----------------------------------------------------------------------*/
double pot(double r, int lang, struct vxaimp *pspa, struct param_pola *ppola)
{
   double v;
   double ionique, coulomb;
   
   vcpot(pspa,r,&ionique,&coulomb);
   v = ionique + coulomb + polapot(ppola,r);
   return(v);
}

