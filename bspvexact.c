

#define NMAX_CBIN 64						//set binomial coefficient array size
double **cbin;

double rgtonorm(double a, int l);
void bbvexact(int nb, int kb, double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx);
//---
double cw3j0(int j1, int j2, int j3);
void set_cbin(double **cbin, int nmax);
double xintegral(int kb, double *tn, int i, int j, int la, int ll, double a, double b);
double xintegral2(int kb, double *tn, double **bspgq, int i, int j, int la, int ll, double a, double b);


/*---------------------------------------------------------------------------
   bbvexact() compute the exchange potential between 2 basis splines, for the
   whole basis set and returns the exchange matrix vvx. The computation is
   done for a given angular momentum value l.
   The exchange operator is the exact form from RGTO.
---------------------------------------------------------------------------*/
void bbvexact(int nb, int kb, double *tn, int lb, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx)
{
   int la,lamax;
   int ll,llmin,llmax;
   struct vshell *shla;						//A copy of the extraction shell for the chosen "l"
   double **rho;						//A copy of density matrix
   int ia,ja,na;
   int ib,jb;
   double wa,suma;
//---
   double **bspgq;
   double t,tmin,tmax,dt;
   int q;

   cbin = dmatrix(0,NMAX_CBIN,0,NMAX_CBIN);
   set_cbin(cbin,NMAX_CBIN);
   bspgq = dmatrix(1,nb,0,QMAX_GAUSS);
   for(ib=1;ib<=nb;ib++)					//set bspline values for Gauss-Legendre integration
      {
      tmin = tn[ib];
      tmax = tn[ib+1+kb];
      dt = tmax - tmin;
      for(q=0;q<QMAX_GAUSS;q++)
        {
        t = (tmin + tmax + dt*egl[q])/2.0;
        bspgq[ib][q] = bspline(tn,bss[ib].ib,kb,t);
        }
      }

   lamax = pspa->lbmax;
   printf("Exchange:\t lamax = %d \n",lamax);

   for(ib=1;ib<=nb;ib++) for(jb=ib;jb<=nb;jb++)
      {
      vvx[ib][jb] = 0.0;
//-------------------------------------
      for(la=0;la<=lamax;la++)
         {
         na = pspa->nlb[la];					//for bookkeeping
         shla = pspa->shlb[la];
         rho = pspa->rho[la];
         llmin = abs(lb-la);
         llmax = lb+la;
         for(ll=llmin;ll<=llmax;ll+=2)
            {
            suma = 0.0;
            for(ia=1;ia<=na;ia++) for(ja=1;ja<=na;ja++)
               {
               wa = rgtonorm(shla[ia].a,la)*rgtonorm(shla[ja].a,la)*rho[ia][ja];
//               suma += wa*xintegral(kb,nbp,tn,ib,jb,la,ll,shla[ia].a,shla[ja].a);
               suma += wa*xintegral2(kb,tn,bspgq,ib,jb,la,ll,shla[ia].a,shla[ja].a);
               }
            vvx[ib][jb] -= (2*la+1)*cw3j0(lb,la,ll)*suma;
//printf("ib=%d \t jb=%d \t lb=%d \t la=%d \t ll=%d \t suma=%+le \t vvx=%+le \n",ib,jb,lb,la,ll,suma,vvx[ib][jb]);
            }
         //vvx[ib][jb] /= 2.0;					//divide by 2 for spin condition (+ or -)
//getchar();
         }
      vvx[ib][jb] /= 2.0;
      vvx[jb][ib] = vvx[ib][jb];				//set transposed element
      }

   //prmatxr2(1,10,vvx);
//-------------------------------------
   //free_dmatrix(cbin,0,NMAX_CBIN,0,NMAX_CBIN);
   //free_dmatrix(bspgq,1,nb,0,QMAX_GAUSS);
}



/*---------------------------------------------------------------------------
   cw3j0(j1,j2,j3) eturns the value of the 3j coefficients:

           | j1  j2  j3 |^2
           | 0   0   0  |

   for a set of 3 values (j1,j2,j3)
---------------------------------------------------------------------------*/
double cw3j0(int j1, int j2, int j3)
{
   double w;
   int g,g2;
   int p1,p2,p3;

   g2 = j1 + j2 + j3;
   if(g2%2 == 0)                                        //even case != 0
      {
      g = g2/2;
      p1 = g-j1;
      p2 = g-j2;
      p3 = g-j3;
      w = cbin[2*p1][p1]*cbin[2*p2][p2]*cbin[2*p3][p3]/((g+1)*cbin[g2+1][g]);
      }
   else
      {
      w = 0.0;
      }

   return(w);
}


/*----------------------------------------------------------------
   set_cbin() returns the binomial coefficients C_nk for n,k <= nmax
----------------------------------------------------------------*/
void set_cbin(double **cbin, int nmax)
{
   int k,n;

   for(n=0;n<=nmax;n++) cbin[n][0] = cbin[n][n] = 1.0;			//get binomial coefficient up to  nmax
   for(n=1;n<=nmax;n++) for(k=1;k<n;k++) cbin[n][k] = cbin[n-1][k] + cbin[n-1][k-1];

//   for(n=0;n<=5;n++) for(k=0;k<=n;k++) printf("n=%d \t k=%d \t cnk=%+le \n",n,k,cbin[n][k]);
//   getchar();
}




double xintegral(int kb, double *tn, int i, int j, int la, int ll, double a, double b)
{
   double vxl;
   double biga,bjgb;
   double t,tmin,tmax,dt;
   double s,smin,smax,ds;
   double sump,sumq;
   int q,p;

   tmin = tn[i];
   tmax = tn[i+1+kb];
   dt = tmax - tmin;
   smin = tn[j];
   smax = tn[j+1+kb];
   ds = smax - smin;

   if(a*tmin*tmin + b*smin*smin - (la+1)*log(tmax*smax) > 36.0) vxl = 0.0;	//integration is useless
   else
      {
      sumq = 0.0;
      for(q=0;q<QMAX_GAUSS;q++)
        {
        t = (tmin + tmax + dt*egl[q])/2.0;
        biga = bspline(tn,i,kb,t)*powint(t,la+1)*exp(-a*t*t);
        sump = 0.0;
        for(p=0;p<QMAX_GAUSS;p++)
           {
           s = (smin + smax + ds*egl[p])/2.0;
           bjgb = bspline(tn,j,kb,s)*powint(s,la+1)*exp(-b*s*s);
           if(s>t) sump += wgl[p]*bjgb*powint(t/s,ll)/s;
	   else sump += wgl[p]*bjgb*powint(s/t,ll)/t;
           }
        sumq += wgl[q]*biga*sump;
        }
      vxl = sumq*dt*ds/4.0;
      }

   return(vxl);
}



double xintegral2(int kb, double *tn, double **bspgq, int i, int j, int la, int ll, double a, double b)
{
   double vxl;
   double biga,bjgb;
   double t,tmin,tmax,dt;
   double s,smin,smax,ds;
   double sump,sumq;
   int q,p;

   tmin = tn[i];
   tmax = tn[i+1+kb];
   dt = tmax - tmin;
   smin = tn[j];
   smax = tn[j+1+kb];
   ds = smax - smin;

   if(a*tmin*tmin + b*smin*smin - (la+1)*log(tmax*smax) > 36.0) vxl = 0.0;	//integration is useless
   else
      {
      sumq = 0.0;
      for(q=0;q<QMAX_GAUSS;q++)
        {
        t = (tmin + tmax + dt*egl[q])/2.0;
        biga = bspgq[i][q]*powint(t,la+1)*exp(-a*t*t);
        sump = 0.0;
        for(p=0;p<QMAX_GAUSS;p++)
           {
           s = (smin + smax + ds*egl[p])/2.0;
           bjgb = bspgq[j][p]*powint(s,la+1)*exp(-b*s*s);
           if(s>t) sump += wgl[p]*bjgb*powint(t/s,ll)/s;
	   else sump += wgl[p]*bjgb*powint(s/t,ll)/t;
           }
        sumq += wgl[q]*biga*sump;
        }
      vxl = sumq*dt*ds/4.0;
      }

   return(vxl);
}
