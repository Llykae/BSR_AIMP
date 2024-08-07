/*-------------------------------------------------------------------------------
------------------------------------------------------------------------------*/

#define LBASMAX 5				//largest angular momentum for basis + 1
#define LEXTMAX 7				//largest angular momentum for extraction + 1

void readorbat(int *nshb, struct shell_basis **shb, int *lbmax, double *zat0, double *qat0, double ***rho0, char *flname);
void initxshell(int *lxmax, int *nlshx, struct shell_basis **shx, char *flname);
void extraction(int lxmax, int *nlshx, struct shell_basis **shx, int lbmax, int *nlshb, struct shell_basis **shbss, double ***rho0);

double rgtonorm(double a, int l);
double rgoverlap(double a, double b, int l);
double mkqlnint(int l, int n, double x, double y);
double mkqlnint1(int ll, int n, double x, double y);
double mkq000int(double x, double y);
double mkq011int(double x, double y);

double cw3j0(int j1, int j2, int j3);



//get orbitals and associated basis shell
void readorbat(int *nlshb, struct shell_basis **shb, int *lbmax, double *zat0, double *qat0, double ***rho0, char *flname)
{
   FILE *fdin;
   char strtmp[1024];
   char *sf;
   double zat,qat;
   double atmp,ctmp;
   double **rho;
   double **corb[LBASMAX];
//   int nbl[10];						//at most l=9
   int nnl[LBASMAX];						//at most l=LBASMAX-1
   int l,lmax;
   int i,j,imax;
   int n,nmax;
   int ish;
   double ai,aj;
   double sum;
//---
   int iprint = 1;

   if(iprint) printf("Orbital basis file name : '%s' \n",flname);

   fdin=fopen(flname,"r");
   for(i=1;i<=4;i++) fgets(strtmp,511,fdin);
   sf = fgets(strtmp,511,fdin);
   sscanf(strtmp,"%le  %d",&zat,&lmax);
   if(iprint) printf("\tZat=%+1.3lf \t lmax=%d \n",zat,lmax);
   *lbmax = lmax;
   *zat0 = zat;
   ish = 0;
   for(l=0;l<=lmax;l++)
      {
      sf = fgets(strtmp,1023,fdin);				//skip 1 title line
//printf("strtmp 1 : %s \n",strtmp);
      sf = fgets(strtmp,1023,fdin);
//printf("strtmp 2 : %s \n",strtmp);
      sscanf(strtmp,"%d %d",&imax,&nmax);
      nlshb[l] = imax;
      nnl[l] = nmax;
      shb[l] = malloc((imax+1)*sizeof(struct shell_basis));	//memory for basis shells for 'l'
printf("l=%d \t imax=%d \t nmax=%d \n",l,imax,nmax);
      sf = fgets(strtmp,1023,fdin);
//printf("strtmp 3 : %s \n",strtmp);
//getchar();
      ish = 0;
      for(i=1;i<=imax;i++)  sf = fgets(strtmp,1023,fdin);		//skip coefficient array
//printf("strtmp  : %s \n",strtmp);
//getchar();}
      }
   fclose(fdin);

//---------------------------------
   for(l=0;l<=lmax;l++)					//set memory for orbitals
      {
      corb[l] = dmatrix(1,nlshb[l],1,nnl[l]);
      }

   fdin=fopen(flname,"r");
   for(i=1;i<=5;i++) sf = fgets(strtmp,255,fdin);
   qat = 0.0;
   for(l=0;l<=lmax;l++)
      {
      sf = fgets(strtmp,255,fdin);
      sf = fgets(strtmp,255,fdin);
      sscanf(strtmp,"%d %d",&imax,&nmax);
//printf("strtmp : %s \n",strtmp);
//getchar();
//printf("l=%d \t nl=%d \n",l,nl);
      ish = 0;
      for(i=1;i<=imax;i++)
         {
         fscanf(fdin,"%lf\n",&atmp);
         ish++;
         shb[l][ish].l = l;
         shb[l][ish].a = atmp;
         shb[l][ish].iat = 1;
	 }
      for(i=1;i<=imax;i++) for(n=1;n<=nmax;n++) fscanf(fdin,"%lf\n",&corb[l][i][n]);
      qat += 2.0*(2*l+1)*nmax;				//cumulate core lectronic charge
      }
   fclose(fdin);
   *qat0 = qat;

   if(iprint) for(l=0;l<=lmax;l++) for(ish=1;ish<=nlshb[l];ish++)
      {
      printf("iat=%d \t l=%d \t ish=%d \t a=%e \t norbl=%d \n",shb[l][ish].iat,shb[l][ish].l,ish,shb[l][ish].a,nnl[shb[l][ish].l]);  
      }

   if(iprint) for(l=0;l<=lmax;l++)
      {
      printf("l=%d \t nbl=%d \t nnl=%d\n",l,nlshb[l],nnl[l]);
      for(i=1;i<=nlshb[l];i++)
         {
         for(n=1;n<=nnl[l];n++) printf("%+le\t",corb[l][i][n]);
         printf("\n");
         }
      printf("\n");
      }

   for(l=0;l<=lmax;l++) 						//loop over all orbital with momentum 'l'
      {
      printf("l=%d \t nbl=%d \t nnl=%d\n",l,nlshb[l],nnl[l]);
      imax = nlshb[l];							//for bookkeeping
      rho0[l] = dmatrix(1,imax,1,imax);					//memory for radial density matrix
      for(n=1;n<=nnl[l];n++)
         {
         sum = 0.0;
         for(i=1;i<=imax;i++) for(j=1;j<=imax;j++)
            {
            ai = shb[l][i].a;
            aj = shb[l][j].a;
            sum += corb[l][i][n]*corb[l][j][n]*rgoverlap(ai,aj,l);
            }
         for(i=1;i<=imax;i++) corb[l][i][n] /= sqrt(sum);		//orbital normalization
//---
         if(iprint)
            { 
            sum = 0.0;
            for(i=1;i<=imax;i++) for(j=1;j<=imax;j++)
               {
               ai = shb[l][i].a;
               aj = shb[l][j].a;
               sum += corb[l][i][n]*corb[l][j][n]*rgoverlap(ai,aj,l);
               }
            printf("\tn=%d \t sum2=%+le \n",n,sum);
            }
         }
      for(i=1;i<=imax;i++) for(j=i;j<=imax;j++)
         {
         rho0[l][i][j] = 0.0;						//matrix density initialization
         for(n=1;n<=nnl[l];n++) rho0[l][i][j] += 2.0*corb[l][i][n]*corb[l][j][n];
         rho0[l][j][i] = rho0[l][i][j];					//forces symmetry
         }
//---
      if(iprint)
         { 
         printf("density matrix from orbitals:\n");
         prmatxr(imax,rho0[l]);
         }
      }

//---------------------------------


//---------------------------------
   for(l=0;l<=lmax;l++)							//free memory for orbitals
      {
      free_dmatrix(corb[l],1,nlshb[l],1,nnl[l]);
      free(shb[l]);
      }
}



/*-----------------------------------------------------------------------------
   Allocate memory for structure shx and reads the basis file to fill the
   structure.
-----------------------------------------------------------------------------*/
void initxshell(int *lxmax, int *nlshx, struct shell_basis **shx, char *flname)
{
   FILE *fdin;
   char strtmp[512];
   char *sf;
   int i,itmp;
   int ish;
   int lmax;
   int nl;
   int l;
   double atmp;

   printf("extraction basis file name : '%s' \n",flname);

   fdin=fopen(flname,"r");
   for(i=1;i<=9;i++) sf = fgets(strtmp,511,fdin);
   sf = fgets(strtmp,511,fdin);
   sscanf(strtmp,"%d",&lmax);
   *lxmax = lmax;

   ish = 0;
   for(l=0;l<=lmax;l++)
      {
      sf = fgets(strtmp,511,fdin);
//printf("A: strtmp : %s \n",strtmp);
      sscanf(strtmp,"%d %d",&itmp,&nl);
//printf("l=%d \t nl=%d \n",l,nl);
      shx[l] = malloc((nl+1)*sizeof(struct shell_basis));		//memory for extraction shells
      nlshx[l] = nl;
      printf("extraction number of radial shells: l=%d \t nlshx=%d \n",l,nlshx[l]);
//      sf = fgets(strtmp,255,fdin);
//printf("B: strtmp : %s \n",strtmp);
      ish = 0;
      for(i=1;i<=nl;i++)
         {
         fscanf(fdin,"%lf\n",&atmp);
         ish++;
         shx[l][ish].l = l;
         shx[l][ish].a = atmp;
         shx[l][ish].iat = 1;
printf("iat=%d \t l=%d \t ish=%d \t a=%e \n",shx[l][ish].iat,shx[l][ish].l,ish,shx[l][ish].a);  
	 }
      }
   fclose(fdin);
   //free(sf);
/*---
   fdout=fopen("xbasis.dat","w");
   for(ibs=1;ibs<=ibsmax;ibs++)
      {
      fprintf(fdout,"%d \t %e \t %e \t %e \t %e \t %e \t %d \t %d \t %d \t %d \t %e \n",ibs,(*bss)[ibs].r[0],(*bss)[ibs].r[1],(*bss)[ibs].r[2], (*bss)[ibs].r[3],(*bss)[ibs].a,(*bss)[ibs].p0,(*bss)[ibs].p1,(*bss)[ibs].p2,(*bss)[ibs].p3,(*bss)[ibs].nrm);  
      }
   fclose(fdout);
---*/
}



/*-------------------------------------------------------------------------------
   extraction() computes the extraction of the exchange potential in the radial
   extraction basis made of GTO. For a given set (a,l) of GTO exponent and angular
   momentum, a single evaluation is done for any magentic number ml.
   Note that the density matrix is normalized to a double ocupation of each
   spatial orbital, while the exchange couples only orbitals with the same spin.

   The extraction is done for each value of lambda separately, as the different
   values are not coupled for a spherically symmetric atom, with closed shell.
------------------------------------------------------------------------------*/
void extraction(int lxmax, int *nlshx, struct shell_basis **shx, int lbmax, int *nlshb, struct shell_basis **shbss, double ***rho0)
{
   double x,y;
   double qxln;
   double qx00;
   double rabcd;
   double clll;
   int l,n;
   int ll,llmin,llmax;
   int lb,lib,ljb;
   int ib,jb;
   int ax,bx,px,qx,pxmax;
//---
   int lambda;
   double **aax,**ssx,**qqx;
   double **kkx,**ssxm1,*ds;
//---
   FILE *fout;

   printf("lbmax=%d \t lxmax=%d \n",lbmax,lxmax);

//-------------------
   fout = fopen("output/vaimp.out","a");
   fprintf(fout,"VEXCHANGE\n");
   fprintf(fout,"# largest angular momentum for extraction\n");
   fprintf(fout,"%d\n",lxmax);
   for(lambda=0;lambda<=lxmax;lambda++)						//extraction in each lambda space
      {
      pxmax = nlshx[lambda];
      printf("# ---\t lambda=%d \t pxmax=%d \n\t",lambda,pxmax);
      aax = dmatrix(1,pxmax,1,pxmax);
      kkx = dmatrix(1,pxmax,1,pxmax);
      for(px=1;px<=pxmax;px++) for(qx=1;qx<=pxmax;qx++)				//loop over extraction shell
         {
         aax[px][qx] = 0.0;
         for(lb=0;lb<=lbmax;lb++) for(ib=1;ib<=nlshb[lb];ib++) for(jb=1;jb<=nlshb[lb];jb++)	//loop over orbital shells
            {
//            lib = ljb = lb;
            n = lambda + lb + 2;
//            printf("ib=%d \t jb=%d \t lnb=%d \t lmb=%d \t rho=%+le \t n=%d \n",ib,jb,lib,ljb,rho1[ib][jb],n);
            llmin = abs(lambda - lb);
            llmax = lambda + lb;
            rabcd = rgtonorm(shbss[lb][ib].a,lb)*rgtonorm(shbss[lb][jb].a,lb);
            rabcd *= rgtonorm(shx[lambda][px].a,lambda)*rgtonorm(shx[lambda][qx].a,lambda);
//            printf("\t rabcd=%+le \n",rabcd);
            for(ll=llmin;ll<=llmax;ll+=2)
               {
               x = shbss[lb][ib].a + shx[lambda][px].a;
               y = shbss[lb][jb].a + shx[lambda][qx].a;
               clll = (2*lb+1)*cw3j0(lb,ll,lambda);
//               if((lb==0) && (ll==0) &&(lambda==0)) qx00 = (2*ll+1)*rabcd*mkq000int(x,y);
//               if((lb==0) && (ll==1) &&(lambda==1)) qx00 = (2*ll+1)*rabcd*mkq011int(x,y);
//               qx00 = rabcd*mkqlnint1(ll,n,x,y);
//               qxln = (2*ll+1)*rabcd*mkqlnint(ll,n,x,y);
               qxln = rabcd*mkqlnint1(ll,n,x,y);
//               printf("\t ll=%d \t x=%+le \t y=%+le \t qxln=%+le \t qx00=%+le\t ratio=%+le \n",ll,x,y,qxln,qx00,qxln/qx00);
//               printf("lb=%d \t ll=%d \t lambda=%d \t clll=%+le \t cw=%+le \n",lb,ll,lambda,clll,cw3j0(ljb,ll,lambda));
               aax[px][qx] += clll*rho0[lb][ib][jb]*qxln;
               }
            }
         aax[px][qx] /= 2.0;						//divide by to because Tr(rho.S) = 2 per shell
//getchar();
         }
//      printf("density matrix: \n");
//      prmatxr(6,rho1);
//      printf("exchange matrix in extraction basis: \n");
//      prmatxr(pxmax,aax);
//-------------------
      ssx = dmatrix(1,pxmax,1,pxmax);
      qqx = dmatrix(1,pxmax,1,pxmax);
      ssxm1 = dmatrix(1,pxmax,1,pxmax);
      ds = dvector(1,pxmax);
      for(px=1;px<=pxmax;px++) for(qx=1;qx<=pxmax;qx++)
         {
         ssx[px][qx] = rgoverlap(shx[lambda][px].a,shx[lambda][qx].a,shx[lambda][px].l);
         }
//      printf("overlap matrix: \n");
//      prmatxr(pxmax,ssx);
      rdiag(pxmax,ssx,ds,qqx);
//      for(px=1;px<=pxmax;px++) printf("px=%d \t ds=%+le \n",px,ds[px]);
      for(px=1;px<=pxmax;px++) for(qx=1;qx<=pxmax;qx++)
         {
         ssxm1[px][qx] = 0.0;
         for(n=1;n<=pxmax;n++) ssxm1[px][qx] += qqx[px][n]*qqx[qx][n]/ds[n];	//get inverse matrix S^-1
         }
//      printf("inverse overlap matrix: \n");
//      prmatxr(pxmax,ssxm1);
//------------------------
//      for(px=1;px<=pxmax;px++) for(qx=1;qx<=pxmax;qx++)
//         {
//         kkx[px][qx] = 0.0;
//         for(ax=1;ax<=pxmax;ax++) kkx[px][qx] += ssxm1[px][ax]*ssx[ax][qx];
//         }
//printf("check inversion:\n");
//prmatxr(pxmax,kkx);
//getchar();
//------------------------
      for(px=1;px<=pxmax;px++) for(qx=1;qx<=pxmax;qx++)
         {
         kkx[px][qx] = 0.0;
         for(ax=1;ax<=pxmax;ax++) for(bx=1;bx<=pxmax;bx++) kkx[px][qx] += ssxm1[px][ax]*aax[ax][bx]*ssxm1[bx][qx];
         }
      printf("extracted exchange matrix: \n");
      prmatxr(pxmax,kkx);
//      fprintf(fout,"--- lambda=%d\n",lambda);
      fprintf(fout,"# lambda\tnlambda\n");
      fprintf(fout,"%d\t%d\n",lambda,pxmax);
      for(px=1;px<=pxmax;px++) fprintf(fout,"%+1.4le           ",shx[lambda][px].a);
      fprintf(fout,"\n");
      for(px=1;px<=pxmax;px++)
         {
         for(qx=1;qx<=pxmax;qx++) fprintf(fout,"%+1.13le  ",kkx[px][qx]);
         fprintf(fout,"\n");
         }


//-------------------
      free_dmatrix(aax,1,pxmax,1,pxmax);
      free_dmatrix(kkx,1,pxmax,1,pxmax);
      free_dmatrix(ssx,1,pxmax,1,pxmax);
      free_dmatrix(qqx,1,pxmax,1,pxmax);
      free_dmatrix(ssxm1,1,pxmax,1,pxmax);
      free_dvector(ds,1,pxmax);
      }

//-------------------
   fclose(fout);
}



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
   rgoverlap() returns the overlap of 2 radial GTO function with exponent 'a'
   and 'b' and total angular momentum 'l' for both.
   NB: the norm of the radial GTO is included.
---------------------------------------------------------------------------*/
double rgoverlap(double a, double b, int l)
{
   double x,s;
   int i;

   x = sqrt(4.0*a*b)/(a + b);
   s = sqrt(x);
   for(i=0;i<=l;i++) s *= x;
   return(s);
}




/*---------------------------------------------------------------------------
   mkqlnint() returns the momentum integration of the gaussian form factor
   product for exchange integrals.
   The product of the norm of the 4 GTO is given as rabcd = Na.Nb.Nc.Nd,
   which was removed from the Iln functions here.

	qln(x,y) = 1/pi Sum dq Iln(x,q).Iln(y,q)

   Test:	l,n = 0,0	OK
   		l,n = 1,1	OK
---------------------------------------------------------------------------*/
double mkqlnint(int ll, int n, double x, double y)
{
   double qlnxy;
   double s,t,u,v,z;
   double a;
   double *mlk,*clij;
   double ui,vj;
   double glp;
   int i,j,k;
   int p;

   if((n - ll - 2)%2 != 0) return (0.0);				//odd case == 0.0

   p = (n - ll - 2)/2;
   z = x + y;
   s = x*y;
   t = s/z;
   u = y/z;
   v = x/z;

   a = sqrt(pi/z)/(8.0*s);
   for(i=1;i<=ll;i++) a /= z;						//divide by z^ll
   for(i=1;i<=p;i++) a /= s;						//divide by s^p
   glp = 1.0;
   for(i=1;i<=p;i++) glp *= (i + ll) + 0.5;				//form Gamma(p+ll+3/2)/Gamma(ll+3/2)
   a *= glp*glp;							//form  prefactor

   mlk = dvector(0,p);
   mlk[0] = 1.0;
   for(i=1;i<=p;i++) mlk[i] = mlk[i-1]*(i-1-p)/(i*(ll+i+0.5));

   clij = dvector(0,2*p+ll);
   clij[0] = 1.0;
   for(i=1;i<=2*p+ll;i++) clij[i] = clij[i-1]*(i-0.5);

   qlnxy = 0.0;
   ui = 1.0;
   for(i=0;i<=p;i++)
      {
      vj = 1.0;
      for(j=0;j<=p;j++)
         {
         qlnxy += ui*vj*mlk[i]*mlk[j]*clij[ll+i+j];
         vj *= v;						//cumulate v^j
         }
      ui *= u;							//cumulate u^i
      }
//   double at;
//   at = sqrt(pi*pi*pi/z)/(16.0*s);
//   printf("ll=%d \t p=%d \t a=%+le \t q=%+le \t at=%+le \n",ll,p,a,qlnxy,at);
   qlnxy *= a;

   free_dvector(mlk,0,p);
   free_dvector(clij,0,2*p+ll);

   return(qlnxy);
}

/*---------------------------------------------------------------------------
   mkqlnint0() returns the momentum integration of the gaussian form factor
   product for exchange integrals. The integration is perform directly and
   expressed in term of Hypergeometric functions 2F1. We have:

	qln(x,y) = Gamma(a+b)/(4b (x+y)^(a+b))
		  {(1+1/s) 2F1(1,-p;b+1;1/s) + (1+s) 2F1(1,-p;b+1;s)}

   with:	1-a = -p
   		a = (n - ll)/2			integer
		b = (n + ll + 1)/2		integer + one half
		s = x/y

   Test:	l,n = 0,0	OK
   		l,n = 1,1	OK
---------------------------------------------------------------------------*/
double mkqlnint1(int ll, int n, double x, double y)
{
   double qlnxy;
   double s,t,u,v,z;
   double q0;
   double b1,cpk,fpk,sk;
   int i,k;
   int a2,b2,p;

   if((n - ll)%2 != 0) return (0.0);					//odd case == 0.0

   a2 = n - ll;								//2a
   b2 = n + ll + 1;							//2b
   b1 = 0.5*b2 + 1.0;
   p = (a2 - 2)/2;							//p = a-1
   z = x + y;
   s = x/y;
   u = 1.0 + s;
   v = 1.0 + 1.0/s;

   q0 = sqrt(pi/z)/(2.0*b2);
   for(i=1;i<=n;i++) q0 *= (i-0.5)/z;

   cpk = 1.0;
   sk = s;
   fpk = cpk*(u + v);
   for(k=1;k<=p;k++)
      {
      cpk *= (p-k)/(b1+k);						//c_k+1 = c_k (p-k)/(b1+k)
      sk *= s;								//s^k			
      fpk += cpk*(u*sk + v/sk);
      }

   qlnxy = q0*fpk;

   return(qlnxy);
}



double mkq000int(double x, double y)
{
   double qlnxy;

   qlnxy = sqrt(pi)/(8.0*x*y*sqrt(x+y));
   return(qlnxy);
}


double mkq011int(double x, double y)
{
   double qlnxy;

   qlnxy = sqrt(pi)/(16.0*x*y*pow(x+y,1.5));
   return(qlnxy);
}

double cw3j0(int j1, int j2, int j3)
{
   double w;
   int g,g2;
   int p1,p2,p3;

   g2 = j1 + j2 + j3;
   if(g2%2 == 0)					//even case != 0
      {
      g = g2/2;
      p1 = g-j1;
      p2 = g-j2;
      p3 = g-j3;
      w = bin[2*p1][p1]*bin[2*p2][p2]*bin[2*p3][p3]/((g+1)*bin[g2+1][g]);
      }
   else
      {
      w = 0.0;
      }

   return(w);
}
