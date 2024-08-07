double mkssab0(double ga, double gb, int la);
void rhonorm(int lb, int nlb, struct vshell *shl, double **rhol);

/*---------------------------------------------------------------------------
   mkssab0() returns the overlap of 2 radial gaussian functions Xa and Xb
   defined as:
        Xa(r) = Na r^la exp(-ga.r^2)
        Xb(r) = Nb r^lb exp(-gb.r^2)
   where both functionis are centerd on the same atomic center.
   Only la == lb is evaluated.
---------------------------------------------------------------------------*/
double mkssab0(double ga, double gb, int la)
{
   double za,zb,zab;
   double sab;
   int il;

   sab = ga + gb;
   za = rgtonorm(ga,la);
   zb = rgtonorm(gb,la);
   zab = 0.5*sqrt(pi/sab)*za*zb;
   for(il=0;il<=la;il++) zab *= (il+0.5)/sab;
   sab = zab;

   return(sab);
}


void rhonorm(int lb, int nlb, struct vshell *shl, double **rhol)
{
   int i,j;
   double **ssl;
   double sum,sumi,sumd,sumu;

   ssl = dmatrix(1,nlb,1,nlb);
   //printf("lb=%d \t nlb=%d \n",lb,nlb);
   //for(i=1;i<=nlb;i++) printf("i=%d \t ai=%+le \n",i,shl[i].a);
   for(i=1;i<=nlb;i++) for(j=1;j<=nlb;j++) ssl[i][j] = mkssab0(shl[i].a,shl[j].a,lb);
   sum = 0.0;
   for(i=1;i<=nlb;i++) for(j=1;j<=nlb;j++) sum += ssl[i][j]*rhol[i][j];
   sumu = ceil(sum);
   sumd = floor(sum);
   if(fabs(sumu-sum) > fabs(sumd-sum)) sumi = sumd;
   else sumi = sumu;
//prmatxr2(1,nlb,1,nlb,rhol);
//prmatxr2(1,nlb,1,nlb,ssl);
//printf("A:\t sum=%+1.16le \t sumd=%+1.16le \t sumu=%+1.16le \t sumi=%+1.16le\n",sum,sumd,sumu,sumi);
   for(i=1;i<=nlb;i++) for(j=1;j<=nlb;j++) rhol[i][j]*=sumi/sum;;
//---
   sum = 0.0;
   for(i=1;i<=nlb;i++) for(j=1;j<=nlb;j++) sum += ssl[i][j]*rhol[i][j];
   sumu = ceil(sum);
   sumd = floor(sum);
   if(fabs(sumu-sum) > fabs(sumd-sum)) sumi = sumd;
   else sumi = sumu;
//printf("B:\t sum=%+1.16le \t sumd=%+1.16le \t sumu=%+1.16le \t sumi=%+1.16le\n",sum,sumd,sumu,sumi);
//getchar();
//---
   free_dmatrix(ssl,1,nlb,1,nlb);
//   exit(1);
}
