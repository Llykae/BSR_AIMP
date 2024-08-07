/*---------------------------------------------------------------------------
  This file contains mathematical function for cartesian Gaussian functions
  generation and some integral computation.

  gto(a,i,j,k,r) = x^i y^j z^k exp(-s) 		with s^2 = a.r^2

  s = (u,v,w) = (bx,by,bz)			with b = sqrt(a);

  B. GERVAIS --- 10/2014
  ---------------------------------------------------------------------------*/
#define pi acos(-1.)
#include "gto.h"



/*---------------------------------------------------------------------------
  gtonorm() returns the norm of the HGTO function for a given exponent a,
  and angular momenta px,py,pz.
  ---------------------------------------------------------------------------*/
double gtonorm(double a, int px, int py, int pz)
{
  double b,nrm;
  int ix,iy,iz;

  b = sqrt(a);
  nrm = pow(2.0*a/pi,0.75);
  for(ix=1;ix<=px;ix++) nrm *= 2.0*b/sqrt(2.0*ix-1.0);
  for(iy=1;iy<=py;iy++) nrm *= 2.0*b/sqrt(2.0*iy-1.0);
  for(iz=1;iz<=pz;iz++) nrm *= 2.0*b/sqrt(2.0*iz-1.0);
  return(nrm);
}


double fnorm(double a, int p1, int p2, int p3)
{
  double b,nrm;
  int i1,i2,i3;

  b = sqrt(a);
  nrm = pow(2*a/pi,0.75);
  for(i1=1;i1<=p1;i1++) nrm *= 2*b/sqrt(2*i1-1.0);
  for(i2=1;i2<=p2;i2++) nrm *= 2*b/sqrt(2*i2-1.0);
  for(i3=1;i3<=p3;i3++) nrm *= 2*b/sqrt(2*i3-1.0);
  return(nrm);
}

/*---------------------------------------------------------------------------
  gtoovrlp() returns the overlap matrix of two series of GTO functions,
  bss1 and bss2.
  ---------------------------------------------------------------------------*/
double gtoovrlp(struct orbital *bss, int i, int j, int k, int h)
{
  double g,fx,fy,fz,c0,pig;

  g = 1.0/(bss[k].arg[i] + bss[h].arg[j]);
  //g2 = 0.25/bss[k].arg[i] + 0.25/bss[h].arg[j];
  pig = pi*g;
  c0 = bss[k].nrm[i]*bss[h].nrm[j]*pig*sqrt(pig);

  fx = g1ovrlp(bss[k].arg[i],bss[h].arg[j],0,0,0.0,0.0);
  fy = g1ovrlp(bss[k].arg[i],bss[h].arg[j],0,0,0.0,0.0);
  fz = g1ovrlp(bss[k].arg[i],bss[h].arg[j],bss[k].l,bss[h].l,0.0,0.0);

  return(c0*fx*fy*fz);
}

/*--------------------------------------------------------------------------
  g1ovrlp: compute the overlap between two gaussian function
  input: exponant mu, nu
  polynomial order m, n;
  position of the gaussian function: xa, xb
  return: overlap
  --------------------------------------------------------------------------*/
double g1ovrlp(double mu, double nu, int m, int n, double xa, double xb)
{
  double eta,x2eta;
  double xp, xpa, xpb;
  double *cp;
  double ss,ss0;
  int pmax, p;
  int ix;

  pmax = m + n;                                                                //order of the polynome on xp
  eta = mu + nu;                                                               //exponant of the resulting gaussian
  x2eta = 2.0*eta;
  xp = (mu*xa + nu*xb)/eta;                                                    //position of the center P
  xpa = xp - xa;                                                               //distance between A and P
  xpb = xp - xb;                                                               //distance between B and P

  cp = dvector(0,pmax);                                                        //coefficient of the polynomial P(xp), memory allocation
  ovrlp_plnmlcoef(m,n,xpa,xpb,cp);                                             //generate the cp coefficient

  ss = 0.0;                                                                    //set to 0.0 the overlap element
  for(p=0;p<=pmax;p+=2)                                                        //loop on even p term only
    {
      ss0 = 1.0;
      for(ix=2;ix<=p;ix+=2) ss0 *= (ix-1.0)/(x2eta);
      ss0 *= cp[p];
      ss += ss0;
    }

  free_dvector(cp,0,pmax);                                                     //free memory
  return(ss);
}



/*--------------------------------------------------------------------
  ovrlp_plnmlcoef: generate the series coefficient cp only for pair
  p values from the expansion  of the product:
  (x+b)^p . (x+c)^q = sum_k a_k x^k
  --------------------------------------------------------------------*/
void ovrlp_plnmlcoef(int m, int n, double xpa, double xpb, double *cp)
{
  double ci[BSLMAX];
  double cj[BSLMAX];
  int p,pmin,pmax;
  int i,j;
  int iswa,iswb;

  if(fabs(xpa)<1e-15) iswa = 0;                                                //set to zero term under threshold
  else iswa = 1;
  if(fabs(xpb)<1e-15) iswb = 0;                                                //set to zero term under threshold
  else iswb = 1;

  pmax = m + n;
  pmin = 0;

  switch(2*iswb + iswa)
    {  
    case 0: 
      pmin = pmax;
      for(p=0;p<pmin;p++) cp[p] = 0.0;                               //loop over p to set all cp[p] identical to 0
      cp[pmax] = 1.0;
      break;

    case 1: 
      pmin = n;
      for(p=0;p<=pmin;p++) cp[p] = 0.0;
      for(i=0;i<=m;i++)
	{  
	  p = pmin + i;
	  cp[p] = bin[m][i]*powint(xpa,m-i);                          //could be faster ...
	}
      break;

    case 2: 
      pmin = m;
      for(p=0;p<=pmin;p++) cp[p] = 0.0;
      for(j=0;j<=n;j++)
	{  
	  p = pmin + j;
	  cp[p] = bin[n][j]*powint(xpb,n-j);
	}
      break;

    case 3:
      for(i=0;i<=m;i++) ci[i] = bin[m][i]*powint(xpa,m-i);
      for(j=0;j<=n;j++) cj[j] = bin[n][j]*powint(xpb,n-j);
      for(p=0;p<=pmax;p++) cp[p] = 0.0;
      pmin = 0;
      for(i=0;i<=m;i++) for(j=0;j<=n;j++)
			  {  
			    p = i + j;
			    cp[p] += ci[i]*cj[j];
			  }                                                                                                                                                                                          
      break;
    default:
      printf("problem in ovrlp_plnmlcoef(): 2iswb + iswa = %d \n",2*iswb + iswa);
      exit(1);
      break;
    }
}


