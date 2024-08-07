double rgtonorm(double a, int l);
double bspgto(int kb,  double *tn, int i, double a, int l);

void vcpot(struct vxaimp *pspa, double r, double *ionique, double *coulomb);
double uclfunc(int l, double q, double r2);
/*---------------------------------------------------------------------------
   rgtonorm() returns the norm of radial GTO function for a given exponent a,
   and total angular momentum l, with l >= -1.
---------------------------------------------------------------------------*/
double rgtonorm(double a, int l)
{
   double b,nrm;
   int i;

   b = 2.0*a;
   nrm = sqrt(4.0*b/pi);
   for(i=0;i<=l;i++) nrm *= b/(i+0.5);
   nrm = sqrt(nrm);
   return(nrm);
}

/*---------------------------------------------------------------------------
   bspgto() returns the overlap of a bspline function Bi(r) with a normalized
   radial GTO function Gla(r) = Nla r^l exp(-a.r^2)
   Integration is done by Gauss-Legendre methode in the definition intervalle
   Bi.A we use here bi(r) = Bi(r)/r

	Iial = Nla Sum_rm^rp bi(r) r^l+1 exp(-a.r^2)
---------------------------------------------------------------------------*/
double bspgto(int kb, double *tn, int i, double a, int l)
{
   double t,tmin,tmax,dt;
   double biga,sum;
   int q;
   int iprn = 0;

   tmin = tn[i];
   tmax = tn[i+1+kb];
   dt = tmax - tmin;
   sum = 0.0;
   if(a*tmin*tmin - (l+1)*log(tmax) > 36.0) sum = 0.0;		//integration is useless
   else for(q=0;q<QMAX_GAUSS;q++)
      {
      t = (tmin + tmax + dt*egl[q])/2.0;
      biga = bspline(tn,i,kb,t)*powint(t,l+1)*exp(-a*t*t);
      sum += wgl[q]*biga;
//printf("t=%e \t bibj=%e \n",t,bibj);
      }
   sum *= rgtonorm(a,l)*dt/2.0;                               //Gaussian normalization here

   if(iprn)
      {
      printf("\n------------------------\n");
      printf("i=%d \t a=%f\t l=%d tmin=%e \t tmax=%e \t sum=%+le \n",i,a,l,tmin,tmax,sum);
      }

   return(sum);
}


void vcpot(struct vxaimp *pspa, double r, double *ionique, double *coulomb)
{
   int la,lamax;
   int nla;
   double zat;
   struct vshell *shla;						//A copy of basis shell for occupied orbitals
   double **rho;						//A copy of density matrix
   double q,r2,vr,vl;
   int i,j;
//---------------

   zat = pspa->zi;
   lamax = pspa->lbmax;

   *ionique = *coulomb = 0.0;

   r2 = r*r;
   *ionique = -zat/r;							//ionic potential
   for(la=0;la<=lamax;la++)					//loop over occupied angular momenta
   {
      nla = pspa->nlb[la];					//for bookkeeping
      shla = pspa->shlb[la];
      rho = pspa->rho[la];
      vl = 0.0;
      for(i=1;i<=nla;i++)					//diagonal element
      {
         q = shla[i].a + shla[i].a;
         vl += rho[i][i]*shla[i].norm*shla[i].norm*uclfunc(la,q,r2);
      }
      vl *= 0.5;
      for(i=1;i<=nla;i++) for(j=i+1;j<=nla;j++)			//off-diagonal elements
      {
         q = shla[i].a + shla[j].a;
         vl += rho[i][j]*shla[i].norm*shla[j].norm*uclfunc(la,q,r2);	//factor 2 for rho(i,j) = rho(j,i)
      }
      *coulomb += (2*la+1)*vl;
   }
}



double uclfunc(int l, double q, double r2)
{
   double x,y,z,ux;

   x = q*r2;
   y = sqrt(x);
   z = erf(y)*sqrt(pi)/(2.0*y);
   switch(l)
      {
      case 0:
	 ux = z/q;
         break;

      case 1:
	 ux = 3.0*z - exp(-x);
         ux /= 2.0*q*q;
         break;

      case 2:
	 ux = 15.0*z - exp(-x)*(2.0*x + 7.0);
         ux /= 4.0*q*q*q;
         break;

      case 3:
	 ux = 105.0*z - exp(-x)*(4.0*x*x + 22.0*x + 57.0);
         ux /= 8.0*q*q*q*q;
         break;

      default:
         printf("no such case yet in uclfunc(): l=%d \n",l);
         exit(1);
      }
   return(ux);
}