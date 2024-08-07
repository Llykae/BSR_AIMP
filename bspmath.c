/*-------------------------------------------------------------------------
   mathematical library for programm lcaomain()
			B. Gervais 01/04
-------------------------------------------------------------------------*/

#define PCSWITCH 0		/*switch for PC compilers, that do not know the functions erf() and erfc()*/




double fact(int n);
void initfact();
void initdblfact();
void logfac();
double gamln(double z);
void intbin();
void intxmg();
void intymg();
void intzpl();
void initnrmzz();
double powint(double z, int k);
double ran2(long *idum);
double gasdev(long *idum);
void lubksb(double **a, int n, int *indx, double *b);
void testlu(double **a, int n, int *indx, double *d);
void prmatxr(int nmax, double **mm);
void prmatxr2(int nmin, int nmax, double **mm);




#if(PCSWITCH)
double erfc(double x);
double erf(double x);
/*---------------------------------------------------------------
   erfc() retruns an estimation of the complementary error function
   with a relative error everywhere less than 1.2e-7
----------------------------------------------------------------*/
double erfc(double x)
{
   double t,z,ans;

   z=fabs(x);
   t=1/(1+0.5*z);
   ans = 0.27886807 + t*(-1.1352039 + t*( 1.48851587 +t*(-0.82215223 + t*0.17087277)));
   ans = -1.26551223 + t*( 1.00002368 + t*( 0.37409196 + t*( 0.09678418 + t*(-0.18628806 + t*ans))));
   ans = t*exp(-z*z + ans);
/*printf("x=%e \t z=%e \t t=%e \t ans=%e \n",x,z,t,ans);*/
   return x>= 0.0 ? ans : 2.0-ans;
}



/*---------------------------------------------------------------
   erf() retruns an estimation of the error function
   with a relative error everywhere less than 1.2e-7
----------------------------------------------------------------*/
double erf(double x)
{
   double ans;
   
   ans = 1.0 - erfc(x);
   return(ans);
}
#endif



/*---------------------------------------------------------------
   erfi() returns an estimation of the complex error function
----------------------------------------------------------------*/
double erfi(double x)
{
   double ans;
   int n;
   double nmax;
   
   nmax = 2*x + sqrt(4*x*x-(4*x-4*40-1));
/*   
printf("x=%e \t nmax=%e \n",x,nmax);
*/
   nmax = 40;
   ans = 0.0;

   for(n=1;n<=nmax;n++)
      {
      ans = ans + 0.5*(exp(-n*n/4.0+n*x)-exp(-n*n/4.0-n*x))/n;
/*
printf("x=%e \t n=%d \t ans=%1.15e \n",x,n,ans);
*/
      }
      
   ans = ans*2.0/pi;
   ans = ans + x/pi; 
/*
printf("x=%e \t ans=%1.15e \n",x,ans);
getchar();
*/ 
   return(ans);
   
}
/*---------------------------------------------------------------
   erfi() returns an estimation of the complex error function
----------------------------------------------------------------*/
double experfi00(double x)
{
   double ans;
   int n;
   double nmax;
   
   nmax = 2*x + sqrt(4*x*x-(4*x-4*40-1));
/*   
printf("x=%e \t nmax=%e \n",x,nmax);
*/
   nmax = 40;
   ans = 0.0;

   for(n=1;n<=nmax;n++)
      {
      ans += (exp(-n*n/4.0+n*x-x*x)-exp(-n*n/4.0-n*x-x*x))/n;
/*
printf("x=%e \t n=%d \t ans=%1.15e \n",x,n,ans);
*/
      }
      
   ans = ans + x*exp(-x*x); 

/*
printf("x=%e \t ans=%1.15e \n",x,ans/pi);
getchar();
*/
   return(ans/pi);
   
}



/*---------------------------------------------------------------
  experfi() returns an estimation of the complex error function erfi(x)
  multiplied by exp(-x*x) and by an additional multiplicative constant
  2/sqrt(pi) i.e. the Dawson's integral.
----------------------------------------------------------------*/
double experfi(double x)
{
  int i,n0;
  double d1,d2,e1,e2,sum,x2,xp,xx;
  double ans;
  static double c1[7];
  static double h;
  static double a1,a2,a3;
  static int init=0;


  if(init == 0)
     {
     init = 1;
     h = 0.4;
     a1 = 2.0/3.0;
     a2 = 0.4;
     a3= 2.0/7.0;
     for(i=1;i<=6;i++)
        {
        c1[i] = exp(-DSQR((2.0*i-1.0)*h));
        }
     }
     
  if(fabs(x) < 0.2)
     {
     x2 = x*x;
     ans = x*(1.0 - a1*x2*(1.0 - a2*x2*(1.0 - a3*x2)));
     }
  else
     {
     xx = fabs(x);
     n0 = 2*(int)(0.5*xx/h+0.5);
     xp = xx - n0*h;
     e1 = exp(2.0*xp*h);
     e2 = e1*e1;
     d1 = n0 + 1;
     d2 = d1 - 2.0;
     sum = 0.0;
     for(i=1;i<=6;i++,d1+=2.0,d2-=2.0,e1*=e2) sum += c1[i]*(e1/d1+1.0/(d2*e1));
     ans = SIGN(exp(-xp*xp),x)*sum/sqpi;
     }

/*     
printf("x=%e \t ans=%1.15e \n",x,2.0*ans/sqpi);
getchar();
*/ 
     

  return(ans*2.0/sqpi);
}


/*---------------------------------------------------------------
  dexperfi() returns an estimation of the complex error function erfi(x)
  multiplied by exp(-x*x) and by an additional multiplicative constant
  2/sqrt(pi) i.e. the Dawson's integral. The precision is better than 1.0 e-15 here.
----------------------------------------------------------------*/
double dexperfi(double x)
{
  int i,n0;
  double d1,d2,e1,e2,sum,x2,xp,xx;
  double ans;
  static double c1[21];
  static double h;
  static double a[11];
  static int init=0;


  if(init == 0)
     {
     init = 1;
     h = 0.3;
     a[0] = 1.0;
     for(i=1;i<=10;i++) a[i] = 1.0/(i + 0.5);
     for(i=1;i<=10;i++) c1[i] = exp(-DSQR((2.0*i-1.0)*h));
     }
     
  if(fabs(x) < 0.2)
     {
     x2 = x*x;
     ans = 1.0;
     for(i=10;i>=1;i--)
        {
        ans = 1.0 - a[i]*x2*ans;
        }
     ans = x*ans;
     }
  else
     {
     xx = fabs(x);
     n0 = 2*(int)(0.5*xx/h+0.5);
     xp = xx - n0*h;
     e1 = exp(2.0*xp*h);
     e2 = e1*e1;
     d1 = n0 + 1;
     d2 = d1 - 2.0;
     sum = 0.0;
     for(i=1;i<=20;i++,d1+=2.0,d2-=2.0,e1*=e2) sum += c1[i]*(e1/d1+1.0/(d2*e1));
     ans = SIGN(exp(-xp*xp),x)*sum/sqpi;
     }

/*     
printf("x=%e \t ans=%1.15e \n",x,2.0*ans/sqpi);
getchar();
*/ 
     

  return(ans*2.0/sqpi);
}



/*---------------------------------------------------------------
  dawson() returns an estimation of the Dawson's integral. The
  precision is better than 1.0 e-15 here.
----------------------------------------------------------------*/
double dawson(double x)
{
   int i,imax,n0;
   double d1,d2,e1,e2,sum,x2,xp,xx;
   double ans;
   static double c1[21];
   static double h;
   static double a[11];
   static double b[15];
   static int init=0;


   if(init == 0)
      {
      init = 1;
      h = 0.3;
      a[0] = 1.0;
      for(i=1;i<=10;i++) a[i] = 1.0/(i + 0.5);
      b[0] = 1.0;
      for(i=1;i<=14;i++) b[i] = i - 0.5;
      for(i=1;i<=20;i++) c1[i] = exp(-DSQR((2.0*i-1.0)*h));
      }
     
   xx = fabs(x);
   if(xx < 0.2)
      {
      x2 = x*x;
      ans = 1.0;
      for(i=10;i>=1;i--)
         {
         ans = 1.0 - a[i]*x2*ans;
         }
      ans = x*ans;
      }

   else if(xx>=0.2 && xx<=10.0)
      {
      n0 = 2*(int)(0.5*xx/h+0.5);
      xp = xx - n0*h;
      e1 = exp(2.0*xp*h);
      e2 = e1*e1;
      d1 = n0 + 1;
      d2 = d1 - 2.0;
      sum = 0.0;
      for(i=1;i<=20;i++,d1+=2.0,d2-=2.0,e1*=e2) sum += c1[i]*(e1/d1+1.0/(d2*e1));
      ans = SIGN(exp(-xp*xp),x)*sum/sqpi;
      }

   else							//xx > 10.0
      {
      x2 = x*x;
      ans = 1.0;
      imax = 12;
      if(xx>33.0) imax = 7;
      if(xx>100.0) imax = 3;
      for(i=imax;i>=1;i--)
         {
         ans = 1.0 + b[i]/x2*ans;
         }
      ans = ans/(2.0*x);
      }
/*     
printf("x=%e \t ans=%1.15e \n",x,2.0*ans/sqpi);
getchar();
*/ 
     

   return(ans);
}



/*---------------------------------------------------------------
  experfi2() returns an estimation of exp(-x*x)*erfi(x) - 1/(x*sqrt(pi))
 
----------------------------------------------------------------*/
double experfi2(double x)
{
  int i,n0;
  double d1,d2,e1,e2,sum,x2,xp,xx;
  double ans;
  static double c1[7];
  static double h;
  static double a1,a2,a3;
  static int init=0;
  double powint();

  if(init == 0)
     {
     init = 1;
     h = 0.4;
     a1 = 2.0/3.0;
     a2 = 0.4;
     a3= 2.0/7.0;
     for(i=1;i<=6;i++)
        {
        c1[i] = exp(-DSQR((2.0*i-1.0)*h));
        }
     }
     
  if(fabs(x) < 0.2)
     {
     x2 = x*x;
     ans = x*(1.0 - a1*x2*(1.0 - a2*x2*(1.0 - a3*x2)));
     ans -= 0.5/x;
     }
  else
     {
     xx = fabs(x);
     n0 = 2*(int)(0.5*xx/h+0.5);
     xp = xx - n0*h;
     e1 = exp(2.0*xp*h);
     e2 = e1*e1;
     d1 = n0 + 1;
     d2 = d1 - 2.0;
     sum = 0.0;
     for(i=1;i<=6;i++,d1+=2.0,d2-=2.0,e1*=e2) sum += c1[i]*(e1/d1+1.0/(d2*e1));
     ans = SIGN(exp(-xp*xp),x)*sum/sqpi;
     ans -= 0.5/x; 
     }

  if(fabs(x) > 10.0)
     {
     ans = 0.25/(powint(x,3)) + 3.0/(8.0*powint(x,5)) + 15.0/(16.0*powint(x,7)) ; 
     }
/*     
printf("x=%e \t ans=%1.15e \n",x,2.0*ans/sqpi);
getchar();
*/ 
     

  return(ans*2.0/sqpi);
}



/*-------------------------------------------------------------------------
   pythag() returns the square root of (a*a + b*b) with an optimal precision.
-------------------------------------------------------------------------*/
double pythag(double a, double b)
{
   double r,fa,fb,x;

   fa = fabs(a);
   fb = fabs(b);

   if(fa > fb)
      {
      x = fb/fa;
      r = fa*sqrt(1.0 + x*x);
      }
   else
      {
      x = fa/fb;
      r = fb*sqrt(1.0 + x*x);
      }
   return(r);
}


/*-------------------------------------------------------------------------
   dot() calculates the scalar product of two vectors a,b whose dimension
   is n.
-------------------------------------------------------------------------*/
double dot(a,b,n)
   double *a,*b;
   int n;
{
   double r=0.;
   int i;

   for(i=1;i<=n;i++) r+=a[i]*b[i];
   return(r);
}



/*-------------------------------------------------------------------------
   prod_mat-vect() calculates the product of matrix a  and vector b whose
   dimension is n.
-------------------------------------------------------------------------*/
void prod_mat_vect(a,b,c,n)

   int n ;
   double *b, *c ;
   double **a ;
 
{
   int i,j ;
  
   for(i=1;i<=n;i++)
      {
      c[i]=0.0 ;
      for(j=1;j<=n;j++)
         {         
         c[i] += a[i][j]*b[j] ;
         }
      }
}


/*-------------------------------------------------------------------------
   prod_mat-vect() calculates the product of matrix a  and vector b whose
   dimension is n.
-------------------------------------------------------------------------*/
void v_mtx_v(double *u, double *v, double **aa, double *s,int n)
{
   int i,j;
   double t;

   *s = 0.0;
   for(i=1;i<=n;i++)
      {
      t = 0.0;
      for(j=1;j<=n;j++)
         {
         t += aa[i][j]*v[j] ;
         }
      (*s) += u[i]*t;
      }
}



/*-------------------------------------------------------------------------
   norm() calculates the  norm of a vector v whose dimension is n.
-------------------------------------------------------------------------*/
void norm(double *v, int n)
{  
   int i;
   double r=0.0;

   for(i=1;i<=n;i++) r+=v[i]*v[i];
   v[0]=sqrt(r);
}


void norm3(double *v)
{  
   double r;

   r = v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
   v[0]=sqrt(r);
}



/*-------------------------------------------------------------------------
   multv calculates the product of a vector b whose dimension is n by
   a constant x
-------------------------------------------------------------------------*/
void multv(x,b,c,n)
   double x,*b,*c;
   int n;
{
   int i;

   for(i=1;i<=n;i++) c[i]=x*b[i];
}



/*----------------------------add()------------------------------
calculates the sum of two vectors a,b whose dimension is n---*/
void vadd(a,b,c,n)
   double *a,*b,*c;
   int n;
{
   int i;

   for(i=1;i<=n;i++) c[i]=a[i]+b[i];
}



/*-------------------------------------------------------------------------
   vdif calculates the difference of two vectors a,b whose dimension is n
-------------------------------------------------------------------------*/
void vdif(double *a, double *b, double *c, int n)
{
   int i;

   for(i=1;i<=n;i++) c[i]=a[i]-b[i];
}


void nvdif(double *a, double *b, double *c, int n)
{
   int i;
   double c2=0.0;

   for(i=1;i<=n;i++)
     {
     c[i]=a[i]-b[i];
     c2 += c[i]*c[i];
     }
   c[0] = sqrt(c2);
}


/*----------------------------ident()------------------------------
     identify two vectors of dimension n
-----------------------------------------------------------------*/
void ident(a,b,n)
   double *a,*b;
   int n;
{
   int i;
   for(i=0;i<=n;i++) b[i]=a[i];
}



/*-------------------------------rot()-------------------------------
    rot calculates rotation from lab to body-centered frame (inv=0)
    and from body centered frame to lab frame (inv=1).call with
    inv=0 is assumed to be first. calculation of inverse only: inv=2

	input  :	r
	output :	rl

      first:  phi= azimuthal angle relative to x axis(positive sense)
      second: thet= polar angle of rl vector relative to x-y plane
                    rotation about intermediate x axis(positive sense)
      third:  psi = azimuthal angle of new x axis
                    relative to intermediate x axis(positive sense)
--------------------------------------------------------------------*/
void rot(rl,psi,thet,phi,inv)
   double *rl,psi,thet,phi;
   int inv;
{
   static double t[4][4];
   double ti[4][4],r[4];
   double a,b,c,d,e,f;
   int i,j;

   ident(rl,r,3);
   switch(inv)
      {
      case 1:	break;

      case 0: 	;

      case 2: 	a = cos(thet);
      		b = sin(thet);
      		c = cos(phi);
      		d = sin(phi);
      		e = cos(psi);
      		f = sin(psi);
      		t[1][1] = a*c*e - d*f;
      		t[1][2] = a*d*e + c*f;
      		t[1][3] = -b*e;
      		t[2][1] = -a*c*f - d*e;
      		t[2][2] = c*e - a*d*f;
      		t[2][3] = b*f;
      		t[3][1] = b*c;
      		t[3][2] = b*d;
      		t[3][3] = a;
		break;
      }
   switch(inv)
      {
      case 0: 	for(i=1;i<4;i++) rl[i]=0;
      		for(i=1;i<4;i++)
		   {
      		   for(j=1;j<4;j++) rl[i]+=t[i][j]*r[j];
		   }
      		break;

      case 1:	;

      case 2: 	for(i=1;i<4;i++)
		   {
		   rl[i]=0;
      		   for(j=1;j<4;j++) ti[j][i]=t[i][j];
		   }
        		for(i=1;i<4;i++)
		   {
      		   for(j=1;j<4;j++) rl[i]+=ti[i][j]*r[j];
		   }
      		break;
      }
}



/*-------------------------------intrpl2()-------------------------------
   two point interpolation of the ordinate y for a given abscissa x
	input:  (xi,yi) abscissa and values of the function to be interpolated
		x value of the function f. Values should be sorted.
	output: y value for y=f(x)
-----------------------------------------------------------------------*/
double intrpl2(x1,x2,z1,z2,x)
   double x1,x2,z1,z2,x;
{
   double r;

   if(x2-x1 == 0.0 ) r=z1;
   else r=z1 + (z2-z1)/(x2-x1)*(x-x1);
   return(r);
}




/*-------------------------------fact()-------------------------------
   returns the value of factorial n! from pre-computed value in array
   fac if n <= NFAC. Stop if the value n exceeds NFAC.
      input: n
      output: f, the factorial of n
-----------------------------------------------------------------------*/
double fact(int n)
{
   if(n <= NFAC) return(fac[n]);
   else
      {
      printf("n=%d too large in fact ... NFAC=%d \n",n,NFAC);
      return(0.0);
      //getchar();
      //getchar();
      }
}



/*------------------------------------------------------------------------
   calculate factorial n for n=0 to NFAC (see globmath)
-----------------------------------------------------------------------*/
void initfact()
{
   int n;
   
   fac[0] = 1;
   for(n=1;n<=NFAC;n++)
      {
      fac[n] = n*fac[n-1];
      }
}



/*------------------------------------------------------------------------
   calculate double factorial n for n=0 to NFAC (see globmath)
   dblfac(n)=(n-1)!!
-----------------------------------------------------------------------*/
void initdblfact()
{
   int i;
   int n;
   int pp;
   
   dblfac[0] = 1.0;

   for(n=1;n<=NFAC;n++)
      {
      if(n%2!=0)
         {
         dblfac[n] = 0.0;
         }
      else
         {
         dblfac[n] = 1.0;
         pp = (n-2)/2;
         for(i=0;i<=pp;i++)
            {
            dblfac[n] *= n-1-2*i;
            }
         }
      }  
}



/*----------------------------logfac()--------------------------
   initialialize a table of the logarithm of the factorial
   called lfac[i] up to NLFAC-1. Both lfac and NLFAC are defined
   in the file globmath which should be included before including
   this file.
----------------------------------------------------------------*/
void logfac()
{
   int i;
   double gamln();

   lfac[0]=0.0;
   for(i=1;i<NLFAC;i++) lfac[i]=gamln(i+1.0);
}



/*--------------------------gamln()----------------------------
   return the logarithm of the gamma function for a real argument.
----------------------------------------------------------------*/
double gamln(double z)
{
   double x,y,t;
   double r=2.5066282746310005,s=1.000000000190015;
   static double cc[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
                -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   int j;

   x=y=z;
   t=x+5.5;
   t-=(x+0.5)*log(t);
   for(j=0;j<=5;j++) s+=cc[j]/++y;
   return(-t+log(r*s/x));
}



/*---------------------------------------------------------------
   intbin() initialize the value of the binomial coefficient for n and k,
   k < n and n,k < NBIN (See globmath).
   	bin[n][k] = n! / k!(n-k)!
----------------------------------------------------------------*/
void intbin()
{
   int n,k;
   for(n=0;n<NBIN;n++) bin[n][0]=bin[n][n]=1.0;
   for(n=1;n<NBIN;n++)
      {
      for(k=1;k<n;k++) bin[n][k]=bin[n-1][k] + bin[n-1][k-1];
      }
}



/*----------------------------intxmg()------------------------------
   Initialize the value of the Xnk coefficient for n and k,
   k < n and n,k < NBIN (See globmath).
   
	xmg[m][k] = (2k)!/((-4.0)^k*(2m-2k)!(2k-m)!);

----------------------------------------------------------------*/
void intxmg()
{
   double fact();
   double powint();
   int m,k,md,kd,kmin,kmax;
   double xk;

   for(m=0;m<NBXMG;m++) 
      {
      md = 2*m;
      kmin = (m+1)/2;
      kmax = m;
      for(k=0;k<kmin;k++) xmg[m][k]=0.0;
      xmg[m][kmin] = fact(2*kmin)/(powint(-4.0,kmin)*fact(md-2*kmin)*fact(2*kmin-m));
      for(k=kmin+1;k<=kmax;k++)
         {
         kd = 2*k;
         xk = -k*(kd-1.0)*(m-k+1)*(md-kd+1.0)/(kd-m-1.0)/(kd-m);
         xmg[m][k] = xk*xmg[m][k-1];

/*
xmg[m][k]=fact(kd)/(powint(-4.0,k)*fact(md-kd)*fact(kd-m));
printf("m=%d \t k=%d \t xmg=%e \t xk=%e \n",m,k,xmg[m][k],xk);
*/

         }
      for(k=kmax+1;k<NBXMG;k++) xmg[m][k]=0.0;
      }
}

/*----------------------------intymg()------------------------------
   Initialize the value of the Ynk coefficient for n and k,
   k < n and n,k < NBYMG (See globmath).

	ymg[m][k] = (2k+1)!/(2(-4.0)^k*(2m-2k)!(2k-m)!);

----------------------------------------------------------------*/
void intymg()
{
   double fact();
   double powint();
   int m,k,md,kmin,kmax;
   double yk;

   for(m=0;m<NBYMG;m++) 
      {
      md = 2*m + 1;
      kmin = (m+1)/2;
      kmax = m;
      for(k=1;k<kmin;k++) ymg[m][k]=0.0;
      ymg[m][kmin] = fact(2*kmin+1)/(2*powint(-4.0,kmin)*fact(md-2*kmin)*fact(2*kmin-m));
      for(k=kmin+1;k<=kmax;k++)
         {
         yk = (k+0.5)/(md-2*k);
         ymg[m][k] = yk*xmg[m][k];
/*
ymg[m][k] = fact(kd)/(2*powint(-4.0,k)*fact(md-2*k)*fact(2*k-m));
printf("m=%d \t k=%d \t ymg=%e \t yk=%e \n",m,k,ymg[m][k],yk);
*/

         }
      for(k=kmax+1;k<NBYMG;k++) ymg[m][k]=0.0;
      }
}


/*----------------------------intzpl()------------------------------
   Initialize the value of the Zpl coefficient for k and m,
   m <= k and m,k < NBYMG (See globmath).
----------------------------------------------------------------*/
void intzpl()
{
   double fact();
   double powint();
   double sgn,sgn1;
   int m,k,mmax;

   for(k=0;k<NBYMG;k++) 
      {
      mmax = k/2;				/*integer division here*/
      if(k%2 == 1) sgn = -1.0;
      else sgn = 1.0;
      for(m=0;m<=mmax;m++)
         {
         if(m%2 == 1) sgn1 *= -1.0;
         else sgn1 = 1.0;
         zpl[k][m] = sgn*exp(lfac[k] - lfac[k-2*m] - lfac[m] - m*log(4.0));
         zpl1[k][m] = sgn1*zpl[k][m];
/*
printf("m=%d \t k=%d \t zpl=%e \t zpl1=%e \t k!=%e \n",m,k,zpl[k][m],zpl1[k][m],(double)(fac[k]));
*/
         }
      for(m=mmax+1;m<NBYMG;m++) zpl[m][k]=0.0;
      }
}



/*--------------------------initnrmzz()----------------------------
   Initialize the value of the Zpl coefficient for k and m,
   m <= k and m,k < NBYMG (See globmath).
----------------------------------------------------------------*/
void initnrmzz()
{
   double a0;
   int l,m;

   for(l=0;l<NFAC;l++) 
      {
      m = 0;
      a0 = (2*l+1)/(2.0*pi);
      nrmzz[l][0] = sqrt(a0/2.0);
/*---
printf("l=%d \t m=%d \t nzz=%e \n",l,m,nrmzz[l][0]);
---*/
      for(m=1;m<=l;m++)
         {
         a0 /= (l-m+1)*(l+m);
         nrmzz[l][m] = sqrt(a0);
/*---
printf("l=%d \t m=%d \t nzz=%e \n",l,m,nrmzz[l][m]);
---*/
         }
      for(m=l+1;m<NFAC;m++) nrmzz[l][m]=0.0;
      }
}



/*----------------------------powint()------------------------------
   Raises a double type number z to an integer power k:
			r = z^k;
----------------------------------------------------------------*/
double powint(double z, int k)
{
   double r=1.0;
   int n;

   if(k<0)
      {
      if(z==0.0) printf("problem in powint: z=%e \n",z);
      z=1.0/z;
      k = abs(k);
      }
   for(n=1;n<=k;n++) r *= z;
   return(r);
}



/*-----------------------------------------------------------------------------
   ran2() returns a uniform deviates between 0.0 and 1.0
-----------------------------------------------------------------------------*/
double ran2(long *idum)
{
   static long iy,ir[98];
   static int iff=0;
   int j;
   int M=714025;
   int IA=1366;
   int IC=150889;

   if(*idum < 0 || iff == 0)
      {
      iff=1;
      if((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
      for(j=1;j<=97;j++)
         {
         *idum=(IA*(*idum)+IC) % M;
         ir[j]=(*idum);
         }
      *idum=(IA*(*idum)+IC) % M;
      iy=(*idum);
      }
   j = 1 + 97.0*iy/M;
   if(j > 97 || j < 1) printf("RAN2: This cannot happen.\n");
   iy=ir[j];
   *idum=(IA*(*idum)+IC) % M;
   ir[j]=(*idum);
   return (double) iy/M;
}



/*------------------------------------------------------------------------------------
   gasdev() returns a normally distributed deviates with zero mean and unit variance.
------------------------------------------------------------------------------------*/
double gasdev(long *idum)
{
   double ran2();
   static int iset=0;
   static double gset;
   double gfac,rsq,v1,v2;

   if(*idum<0) iset=0;
   if(iset==0)
      {
      do
         {
         v1=2.0*ran2(idum)-1.0;
         v2=2.0*ran2(idum)-1.0;
         rsq=v1*v1+v2*v2;
         }
      while (rsq>=1.0 || rsq==0);
      gfac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*gfac;
      iset=1;
      return(v2*gfac);
      }
   else
      {
      iset=0;
      return(gset);
      }
} 



/*------------------------------------------------------------------------------------
   ludcmp() performs the LU decomposition of a given square matrix a. See Numerical
   Recipies in C.

				A = L.U

	input:	a	the matrix to decompose
		n	the rank of the matrix

	output:	a	the decomposed matrix in place of the original matrix
		indx	a vector that records the rows permutations
		*d	+/- according to the parity of row interchanges
------------------------------------------------------------------------------------*/
void ludcmp(double **a, int n, int *indx, double *d)
{
   int i,imax,j,k;
   double big,dum,sum,temp;
   double *vv;
   double tiny=1.0e-20;

   vv=dvector(1,n);
   *d = 1.0;
 
   for(i=1;i<=n;i++)
      {
      big = 0.0;
      for(j=1;j<=n;j++) if((temp=fabs(a[i][j])) > big) big = temp;
      if(big == 0.0)
         {
         printf("singular matrix in ludcmp: big=%e \n",big);
         exit(1);
         }
      vv[i] = 1/big;
      }

   for(j=1;j<=n;j++)
      {
      for(i=1;i<j;i++)
         {
         sum = a[i][j];
         for(k=1;k<i;k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         }
      big = 0.0;
      for(i=j;i<=n;i++)
         {
         sum = a[i][j];
         for(k=1;k<j;k++) sum -= a[i][k]*a[k][j];
         a[i][j] = sum;
         if((dum=vv[i]*fabs(sum)) >= big)
            {
            big = dum;
            imax = i;
            }
         }
      if(j != imax)				/*swap columns imax and j*/
         {
         for(k=1;k<=n;k++)
            {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
            }
         *d = -(*d);
         vv[imax] = vv[j];
         }
      indx[j] = imax;
      if(a[j][j] == 0.0) a[j][j] = tiny;
      if(j != n)
         {
         dum = 1.0/a[j][j];
         for(i=j+1;i<=n;i++) a[i][j] *= dum;
         }
      }

   free_dvector(vv,1,n);
}



/*------------------------------------------------------------------------------------
   lubksb() performs the substituion for solving a linear set of equations after
   LU decomposition of a given square matrix a and a given RHS vector b.
   See Numerical Recipies in C.

				A.x = b

	input:	a	the matrix to decompose
		n	the rank of the matrix

	output:	a	the decomposed matrix in place of the original matrix
		indx	a vector that records the rows permutations
		*d	+/- according to the parity of row interchanges
------------------------------------------------------------------------------------*/
void lubksb(double **a, int n, int *indx, double *b)
{
   int i,ii=0,ip,j;
   double sum;

   for(i=1;i<=n;i++)
      {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if(ii) for(j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
      else if(sum) ii = i;
      b[i] = sum;
      }

   for(i=n;i>=1;i--)
      {
      sum = b[i];
      for(j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
      b[i] = sum/a[i][i];
      }

}


void testlu(double **a, int n, int *indx, double *d)
{
   double **y,**c;
   double *col,tij;
   int i,j,k;

   c=dmatrix(1,n,1,n);
   y=dmatrix(1,n,1,n);
   col=dvector(1,n);


printf("a: in testlu ...");
//getchar();
   for(j=1;j<=n;j++)
      {
      for(i=1;i<=n;i++) c[i][j] = a[i][j];
      }

printf("b: in testlu ...");
//getchar();

   ludcmp(c,n,indx,d);

printf("c: in testlu ...");
//getchar();

   for(j=1;j<=n;j++)
      {
      for(i=1;i<=n;i++) col[i] = 0.0;
      col[j] = 1.0;
//printf("j=%d \n",j);
//getchar();
      lubksb(c,n,indx,col);
      for(i=1;i<=n;i++) y[i][j] = col[i];
      }

   for(i=1;i<=n;i++)
      {
      for(j=1;j<=n;j++)
         {
         tij = 0.0;
         for(k=1;k<=n;k++)
            {
            tij += y[i][k]*a[k][j];
            }
         //printf("i=%d \t j=%d \t y=%e \t a=%e \t tij=%e \n",i,j,y[i][j],a[i][j],tij);
         }
//getchar();
      }
}


void prmatxr(int nmax, double **mm)
{
   double rcut = 1.0e-13;
   double r;
   int p,q;

   printf("--------------- mm part:\n");
   for(p=1;p<=nmax;p++)
      {
      for(q=1;q<=nmax;q++)
         {
         r = mm[p][q];
         if(fabs(r) <= rcut) r = 0.0;
         printf("%+1.3e \t",r);
         }
      printf("\n");
      }
}


void prmatxr2(int nmin, int nmax, double **mm)
{
   double rcut = 1.0e-13;
   double r;
   int p,q;

   printf("--------------- mm part:\n");
   for(p=nmin;p<=nmax;p++)
      {
      for(q=nmin;q<=nmax;q++)
         {
         r = mm[p][q];
         if(fabs(r) <= rcut) r = 0.0;
         printf("%+1.3e \t",r);
         }
      printf("\n");
      }
}

