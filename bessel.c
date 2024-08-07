/*-------------------------------------------------------------------------
   This file contains some routines to calculate bessel functions. All of
   them are adaptated from the programm given in Numerical Recipies in C
   (see chapter 5 and 6). We also add a few routines for the legendre
   polynimials and legendre function of the second kind. 
   The following routines are included:
	bsljy:	bessel function of any order.
	bslsph:	spherical bessel function of order l (up to order 60)
	bslrt:	ratio Ja(x)/Ja-1(x)
        bslchb:	G1 and G2 as the Temme serie for routine bsljy
	chebev:	chebyshev expansion of G1 and G2 for the routine bslchb
	plgndr:	Legendre polynomials
------------------------------------------------------------------------*/

#define EPS 1.e-16
#define FPMIN 1.e-32
#define MAXIT 10000
#define XMIN 2.0
#define NUSE1 7
#define NUSE2 8


/*-------------------------------------------------------------------------
   Calculation of the bessel function Jv and Yv and their derivatives for
   v >= 0.
	input:	x	argument of the functions (strictly positive real)
		v	order of the function (positive real)
	output:	*jv	Jv(x)
		*yv	Yv(x)
		*djv	Jv'(x)
		*dyv	Yv'(x)
------------------------------------------------------------------------*/
void bsljy(v,x,jv,yv,djv,dyv)
   double v,x,*jv,*yv,*djv,*dyv;
{
   void bslchb();
   int i,isign,l,nl;
   double x2,xi,xi2,v1,v12,piv1,piv2;
   double a,b,br,bi,c,cr,ci,d,dr,di,del,del1,dlr,dli,den,e,f,ff,h,t;
   double p,q,r,w,sum,sum1;
   double fact,fact2,fact3,gam,gam1,gam2,gammi,gampl;
   double j1,jp1,jl,jpl,jt;
   double ry1,ryv,rjv,ryvp;

   if(x<=0 || v<0) return;				/*x or v out of range*/
   nl=(x<XMIN ? (int)(v+0.5) : MAX(0,(int)(v-x+1.5)));	/*number of reccurences for Jv and Yv)*/
   v1=v-nl;
   v12=v1*v1;
   xi=1/x;
   xi2=2*xi;
   w=xi2/pi;						/*wronskian*/
   isign=1;
   h=v*xi;						/*evaluation of CF1 by modified Lenz method*/
   if(h<FPMIN) h=FPMIN;
   c=h;
   b=v*xi2;
   d=0;
   for(i=1;i<=MAXIT;i++)				/*Lenz reccurence*/
      {
      b+=xi2;
      d=b-d;
      if(fabs(d) < FPMIN) d=FPMIN;
      d=1/d;
      c=b-1/c;
      if(fabs(c) < FPMIN) c=FPMIN;
      del=c*d;
      h*=del;
      if(d<0) isign=-isign;
      if(fabs(del-1)<EPS) break;
      }
   if(i>MAXIT) return;					/*iteration failed*/
   j1=jl=isign*FPMIN;					/*initialize Jv for dowward reccurence*/
   jp1=jpl=h*jl;					/*initialize Jv' for dowward reccurence*/
   fact=v*xi;
   for(l=nl;l>=1;l--)					/*downward reccuresnce*/
      {
      jt=fact*jl+jpl;
      fact-=xi;
      jpl=fact*jt-jl;
      jl=jt;
      }
   if(jl==0) jl=EPS;
   f=jpl/jl;						/*jv'/Jv*/
   if(x<XMIN)						/*use the serie to normalize*/
      {
      x2=0.5*x;
      piv1=v1*pi;
      fact=(piv1<EPS ? 1. : piv1/sin(piv1));
      d=-log(x2);
      e=d*v1;
      fact2=(fabs(e)<EPS ? 1. : sinh(e)/e);
      bslchb(v1,&gam1,&gam2,&gampl,&gammi);		/*chebyshev evaluation of gam1 and gam2*/
      ff=2/pi*fact*(gam1*cosh(e) + gam2*fact2*d);
      e=exp(e);
      p=e/(pi*gampl);
      q=1/(pi*e*gammi);
      piv2=0.5*piv1;
      fact3=(fabs(piv2)<EPS ? 1. : sin(piv2)/piv2);
      r=pi*piv2*fact3*fact3;
      c=1;
      d=-x2*x2;
      sum=ff+r*q;
      sum1=p;
      for(i=1;i<=MAXIT;i++)
         {
         ff=(i*ff+p+q)/(i*i+v12);
         c*=d/i;
         p/=i-v1;
         q/=i+v1;
         del=c*(ff+r*q);
         sum+=del;
         del1=c*p-i*del;
         sum1+=del1;
         if(fabs(del) < (1+fabs(sum))*EPS) break;
         }
      if(i>MAXIT) return;				/*iteration failed*/
      ryv=-sum;
      ry1=-sum1*xi2;
      ryvp=v1*xi*ryv-ry1;
      rjv=w/(ryvp-f*ryv);
      }
   else							/*evaluate CF2 by modified Lenz method*/
      {
      a=0.25-v12;
      p=-0.5*xi;
      q=1;
      br=2*x;
      bi=2;
      fact=a*xi/(p*p+q*q);
      cr=br+q*fact;
      ci=bi+p*fact;
      den=br*br+bi*bi;
      dr=br/den;
      di=-bi/den;
      dlr=cr*dr-ci*di;
      dli=cr*di+ci*dr;
      t=p*dlr-q*dli;
      q=p*dli+q*dlr;
      p=t;
      for(i=2;i<=MAXIT;i++)
         {
         a+=2*(i-1);
         bi+=2;
         dr=a*dr+br;
         di=a*di+bi;
         if(fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
         fact=a/(cr*cr+ci*ci);
         cr=br+cr*fact;
         ci=bi-ci*fact;
         if(fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
         den=dr*dr+di*di;
         dr/=den;
         di/=-den;
         dlr=cr*dr-ci*di;
         dli=cr*di+ci*dr;
         t=p*dlr-q*dli;
         q=p*dli+q*dlr;
         p=t;
/*printf("i=%d \t p=%e \t q=%e \t w=%e \n",i,p,q,w);*/
         if(fabs(dlr-1)+fabs(dli) < EPS) break;
         }
      if(i>MAXIT) return;				/*iteration failed*/
      gam=(p-f)/q;
      rjv=sqrt(w/((p-f)*gam+q));
      rjv=SIGN(rjv,jl);
      ryv=rjv*gam;
      ryvp=ryv*(p+q/gam);
      ry1=v1*xi*ryv-ryvp;
      }
   fact=rjv/jl;
   *jv=j1*fact;						/*scales original Jv*/
   *djv=jp1*fact;					/*scales original Jv'*/
   for(i=1;i<=nl;i++)					/*upward reccurence for Yv*/
      {
      t=(v1+i)*xi2*ry1-ryv;
      ryv=ry1;
      ry1=t;
      }
   *yv=ryv;						/*value of Yv*/
   *dyv=v*xi*ryv-ry1;					/*value of Yv'*/
}


/*-------------------------------------------------------------------------
   Calculation of the modified bessel function Iv and Kv and their 
   derivatives for v >= 0.
	input:	x	argument of the functions (strictly positive real)
		v	order of the function (positive real)
	output:	*iv	Iv(x)
		*kv	Kv(x)
		*div	Iv'(x)
		*dkv	Kv'(x)
------------------------------------------------------------------------*/
void bslik(v,x,iv,kv,div,dkv)
   double v,x,*iv,*kv,*div,*dkv;
{
   void bslchb();
   int i,l,nl;
   double x2,xi,xi2,v1,v12,piv1;
   double a,a1,b,c,d,del,del1,delh,dels,e,f,ff,h,t;
   double p,q,q1,q2,qn,s,sum,sum1;
   double fact,fact2,gam1,gam2,gammi,gampl;
   double ril,ril1,rip1,ripl;
   double rk1,rkv1,rkv1p,riv1;

   if(x<=0 || v<0) return;				/*x or v out of range*/
   nl=(x<XMIN ? (int)(v+0.5) : MAX(0,(int)(v-x+1.5)));	/*number of reccurences for Jv and Yv)*/
   v1=v-nl;
   v12=v1*v1;
   xi=1/x;
   xi2=2*xi;
   h=v*xi;						/*evaluation of CF1 by modified Lenz method*/
   if(h<FPMIN) h=FPMIN;
   c=h;
   b=v*xi2;
   d=0;
   for(i=1;i<=MAXIT;i++)				/*Lenz reccurence*/
      {
      b+=xi2;
      d=1/(b+d);					/*denominator cannot vanish*/
      c=b+1/c;
      del=c*d;
      h*=del;
      if(fabs(del-1)<EPS) break;
      }
   if(i>MAXIT) return;					/*iteration failed*/
   ril1=ril=FPMIN;					/*initialize Iv for dowward reccurence*/
   rip1=ripl=h*ril;					/*initialize Iv' for dowward reccurence*/
   fact=v*xi;
   for(l=nl;l>=1;l--)					/*downward reccuresnce*/
      {
      t=fact*ril+ripl;
      fact-=xi;
      ripl=fact*t+ril;
      ril=t;
      }
   f=ripl/ril;						/*Iv'/Iv*/
   if(x<XMIN)						/*use the serie to normalize*/
      {
      x2=0.5*x;
      piv1=v1*pi;
      fact=(piv1<EPS ? 1. : piv1/sin(piv1));
      d=-log(x2);
      e=d*v1;
      fact2=(fabs(e)<EPS ? 1. : sinh(e)/e);
      bslchb(v1,&gam1,&gam2,&gampl,&gammi);		/*chebyshev evaluation of gam1 and gam2*/
      sum=ff=fact*(gam1*cosh(e) + gam2*fact2*d);
      e=exp(e);
      p=0.5*e/gampl;
      q=0.5/(e*gammi);
      c=1;
      d=x2*x2;
      sum1=p;
      for(i=1;i<=MAXIT;i++)
         {
         ff=(i*ff+p+q)/(i*i-v12);
         c*=d/i;
         p/=i-v1;
         q/=i+v1;
         del=c*ff;
         sum+=del;
         del1=c*(p-i*ff);
         sum1+=del1;
         if(fabs(del) < fabs(sum)*EPS) break;
         }
      if(i>MAXIT) return;				/*iteration failed*/
      rkv1=sum;
      rk1=sum1*xi2;
      }
   else							/*evaluate CF2 by modified Lenz method*/
      {							/*and evaluate S by Thompson-Barnett method*/
      b=2*(1+x);
      h=delh=d=1/b;
      q1=0;
      q2=1;
      a1=0.25-v12;
      q=c=a1;
      a=-a1;
      s=1+q*delh;
      for(i=2;i<=MAXIT;i++)
         {
         a-=2*(i-1);
         c=-a*c/i;
         qn=(q1-b*q2)/a;
         q1=q2;
         q2=qn;
         q+=c*qn;
         b+=2;
         d=1/(b+a*d);
         delh*=b*d-1;
         h+=delh;
         dels=q*delh;
         s+=dels;
         if(fabs(dels/s) < EPS) break;			/*test only the sum convergence*/
         }
      if(i>MAXIT) return;				/*iteration failed*/
      h=a1*h;
      rkv1=sqrt(pi/(2*x))*exp(-x)/s;			/*rescaling may be obtain by skipping exp(-x)*/
      rk1=rkv1*xi*(v1+x+0.5-h);
      }
   rkv1p=v1*xi*rkv1-rk1;
   riv1=xi/(f*rkv1-rkv1p);				/*get Iv1 from the wronskian*/
   *iv=(riv1*ril1)/ril;					/*scales original Iv*/
   *div=(riv1*rip1)/ril;				/*scales original Iv'*/
   for(i=1;i<=nl;i++)					/*upward reccurence for Kv*/
      {
      t=(v1+i)*xi2*rk1-rkv1;
      rkv1=rk1;
      rk1=t;
      }
   *kv=rkv1;						/*value of Yv*/
   *dkv=v*xi*rkv1-rk1;					/*value of Yv'*/
}

	
/*-------------------------------------------------------------------------
   calculation of the spherical bessel function of the first and
   second kind. This calculation is done by reccurence from a
   generating function f(l,z) using the continued fraction.
   input:       l       order (l<60)
		z       argument (z!=0)
   output:      *jl     spherical bessel function of the first kind
		*nl     spherical bessel function of the second kind
		*djl    derivative of jl
		*dnl    derivative of nl
------------------------------------------------------------------------*/
void bslsph(l,z,jl,nl,djl,dnl)
   double z,*jl,*nl,*djl,*dnl;
   int l;
{
   int i;
   double n[100],f,g,b,c,d,x;
   double eps=1e-12;


   n[0]=-cos(z)/z;
   n[1]=-cos(z)/z/z-sin(z)/z;
   for(i=1;i<=l;i++) n[i+1]=(2*i+1)*n[i]/z - n[i-1];
   if(l==0) g=-tan(z) - 1/z;
   else g=n[l-1]/n[l] - (l+1)/z;
   if(n[l]==0) printf("nl=%e \n",n[l]);
   if(l==0) f=c=eps*eps/z;
   else f=c=l/z;
   d=x=0;
   i=1;
   while(fabs(x-1) > eps)
      {
      b=(2*(l+i)+1)/z;
      d=(b-d);
      if(d==0) d=eps*eps;
      d=1/d;
      c=b-1/c;
      if(c==0) c=eps*eps;
      x=c*d;
      f*=x;
      i++;
      }
   *nl=n[l];
   *dnl=g*n[l];
   if(g==f) printf("g=%e \t f=%e \n",g,f);
   *jl=1/z/z/n[l]/(g-f);
   *djl=f*(*jl);
}


/*-------------------------------------------------------------------------
   calculation of the ratio Ja(x)/Ja-1(x) for a real value x and a complex
   value of a.
	input:	ar,ai	real and imaginary part of a
		x	argument of the Bessel function
	output:	rr,ri	real and imagnary value of the ratio.

  We use a formulae given by Abramowitz (See Abramowitz 9-1-73 page363)
  to compute the continued fraction equal to the ratio. We use the Lenz
  method (see Numerical Recipies in C) adapted for complex analysis.
------------------------------------------------------------------------*/
void bslrt(ar,ai,z,rr,ri)
   double ar,ai,z,*rr,*ri;
{
   double fi,fr,gi,gr,br,bi,cr,ci,cm2,dr,di,dm2,s,x,y;
   double eps=1e-12,tiny=1e-24;
   int i;

   if(z==0)
      {
      *rr=0;
      *ri=0;
      return;
      }
   fi=fr=cr=ci=eps*eps;
   cm2=ci*ci+cr*cr;
   x=y=dr=di=0;
   i=1;
   bi=ai/z;
   while(fabs(sqrt(x*x+y*y)-1) > eps)
      {
      if(i==1) s=-1;
      else s=1;
      br=2*(ar+i-1)/z;
      dr=br-s*dr;
      di=bi-s*di;
      dm2=dr*dr+di*di;
      if(dm2==0) dm2=tiny;
      dr=dr/dm2;
      di=-di/dm2;
      cr=br-s*cr/cm2;
      ci=bi+s*ci/cm2;
      cm2=ci*ci+cr*cr;
      if(cm2==0) cm2=tiny;
      x=dr*cr-di*ci;
      y=dr*ci+di*cr;
      gr=x*fr-y*fi;
      gi=x*fi+y*fr;
      fr=gr;
      fi=gi;
      i++;
      }
   *rr=fr;
   *ri=fi;
}


/*-------------------------------------------------------------------------
   bslchb evaluates g1 and g2 by Chebyshev expansion for |x| <= 0.5. It also
   returns 1/G(1+x) and 1/(G(1-x).
	input:	x	argument
	output:	g1
		g2
		gpl	1/G(1+x)
		gmi	1/G(1-x)
------------------------------------------------------------------------*/
void bslchb(x,g1,g2,gpl,gmi)
   double x,*g1,*g2,*gpl,*gmi;
{
   double chebev();
   double xx;
   static double ch1[]={-1.142022680371172e+0, 6.516511267076e-3, 3.08709017308e-4,-3.470626964e-6, 6.943764e-9,3.6780e-11,-1.36e-13};
   static double ch2[]={ 1.843740587300906e+0,-7.6852840844786e-2, 1.271927136655e-3,-4.971736704e-6,-3.3126120e-8, 2.42310e-10,-1.70e-13,-1.0e-15};

   xx=8*x*x-1;					/*make range between -1 and 1*/
   *g1=chebev(-1.0,1.0,ch1,NUSE1,xx);		/*Chebyshev expansion*/
   *g2=chebev(-1.0,1.0,ch2,NUSE2,xx);		/*Chebyshev expansion*/
   *gpl=*g2-x*(*g1);
   *gmi=*g2+x*(*g1);
}


/*---------------------------------------------------------------------
   caluclation of the chebyshev approximation of a function. The
   coefficients of interpolation are given as input (See numerical
   Recipies in C 5.8 page 193).
	input:	a,b	limit of the interval of interpolation
		c	table of coefficient
		m	size of table c
		x	argument of the function (between a and b)
---------------------------------------------------------------------*/
double chebev(a,b,c,m,x)
   double a,b,*c,x;
   int m;
{
   double d=0,dd=0,sv,y,y2;
   int j;

   y=(2*x-a-b)/(b-a);
   y2=2*y;
   for(j=m-1;j>=1;j--)
      {
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
      }
   y=y*d-dd+0.5*c[0];
   return(y);
}


/*-----------------------------------------------------------------------
   bslj0 calculates j0 for any real x by expansion as polynomials of z=8/x
   (x > 8) or rationnal functions (x < 8).
	input:	x	argument
-----------------------------------------------------------------------*/
float bslj0(x)
   double x;
{
   double ax,z,xx,y,ans,ans1,ans2;

   ax=fabs(x);
   if(ax < 8)
      {
      y=x*x;
      ans1=57568490574.0 + y*(-13362590354.0 + y*(651619640.7 + y*(-11214424.18 + y*(77392.33017 - y*184.9052456))));
      ans2=57568490411.0 + y*( 1029532985.0 + y*(9494680.718 + y*(59272.64853 + y*(267.8532712 + y))));
      ans=ans1/ans2;
      }
   else
      {
      z=8/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0 + y*(-1.098628627e-3 + y*(2.734510407e-5 + y*(-2.073370639e-6 + y*2.093887211e-7)));
      ans2=-1.562499995e-2 + y*(1.430488765e-4 + y*(-6.911147651e-6 + y*(7.621095161e-7 - y*9.34935152e-8)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
      }
   return(ans);
}


/*-----------------------------------------------------------------------
   bslj1 calculates j1 for any real x by expansion as polynomials of z=8/x
   (x > 8) or rationnal functions (x < 8).
	input:	x	argument
-----------------------------------------------------------------------*/
float bslj1(x)
   double x;
{
   double ax,z,xx,y,ans,ans1,ans2;

   ax=fabs(x);
   if(ax < 8)
      {
      y=x*x;
      ans1=x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1 + y*(-2972611.439 + y*(15704.48260 - y*30.16036606)))));
      ans2=144725228442 + y*(2300535178.0 + y*(18583304.74 + y*(99447.43394 + y*(376.9991397 + y))));
      ans=ans1/ans2;
      }
   else
      {
      z=8/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0 + y*(1.83105e-3 + y*(-3.516396496e-5 + y*(2.457520174e-6 - y*2.40337019e-7)));
      ans2=4.687499995e-2 + y*(-2.002690873e-4 + y*(8.449199096e-6 + y*(-8.8228987e-7 + y*1.05787412e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
      }
   if(x<0) ans=-ans;
   return(ans);
}



/*--------------------------------------------------------------------------------------
   plgndr computes the legendre polynomial by an upward reccurence relation.
	input:	n 	polynomial order
		x	polynamial argument part
	output:	pn	polynomial value
--------------------------------------------------------------------------------------*/
double plgndr(n,x)
   double x;
   int n;
{
   int i;
   double a,b,p0,p1,p2;

   p0=1;
   p1=x;
   for(i=1;i<n;i++)                                  /*get legendre polynomials pn and pn-1 by reccurence*/
      {
      a=(2*i+1.)/(i+1.);
      b=i/(i+1.);
      p2=a*x*p1-b*p0;
      p0=p1;                                       /*p0 == pn-1*/
      p1=p2;                                       /*p1 == pn*/
      }
   if(n==0) p1=1;
   return(p1);
}


/*--------------------------------------------------------------------------------------
   lgndr computes the legendre polynomial by an upward reccurence relation.
   It returns both the legendre polynomial and its derivative. See Abramowitz chapter 8.
	input:	n 	polynomial order
		x	polynamial argument part
	output:	pn	polynomial value
--------------------------------------------------------------------------------------*/
void lgndr(n,x,pn,dpn)
   double x,*pn,*dpn;
   int n;
{
   int i,nmax;
   double a,b,p0,p1,p2,u;
   double ak[21];

   p0=1;
   p1=x;
   for(i=1;i<n;i++)				/*get legendre polynomials pn and pn-1 by reccurence*/
      {
      a=(2*i+1.)/(i+1.);
      b=i/(i+1.);
      p2=a*x*p1-b*p0;
      p0=p1;					/*p0 == pn-1*/
      p1=p2;					/*p1 == pn*/
      }
   if(n==0)
      {
      *pn=1;
      *dpn=0;
      }
   else
      {
      *pn=p1;
      if(fabs(x*x-1) > 0.05/n)
         {
         *dpn=n*(x*p1-p0)/(x*x-1); 
         }
      else					/*uses hypergeometric expansion around 1*/
         {
         nmax=MIN(n,20);
         u=0.5*(1-fabs(x));
         ak[0]=-0.5;
         ak[1]=0.5*n*(n+1);
         for(i=1;i<nmax;i++)
            {
            ak[i+1]=(-n+i)*(n+i+1.)/((i+1.)*i)*ak[i];
            }
         p2=ak[nmax];
         for(i=nmax-1;i>=1;i--)
            {
            p2=ak[i] + p2*u;;
            }
         if(x>0) *dpn=p2;
         else
            {
            if(n%2==0) *dpn=-p2;		/*n even, dpn/dx antisymetric*/
            else *dpn=p2;			/*n odd, dpn/dx symetric*/
            }
         }
      }
}
