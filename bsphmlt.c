/*---------------------------------------------------------------------
   This file contains the function for hamiltonian matrix eelements
   evaluation for basis splines.
---------------------------------------------------------------------*/

void mkss(int nbasis, int kb, double *tn, double **ss, struct bsp_basis *bss);
void mkttbloch(int nbasis, int lang, int kb, double *tn, double **tt, struct bsp_basis *bss);
void mkvvl(int nbasis, int kb, double *tn, int lang, struct bsp_basis *bss, struct vxaimp *pspa, struct param_pola *ppola, double **pola, double **ionique, double **coulomb);

void mkvvex(int nb, int kb, double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx);
void mkvvpa(int nb, int kb,  double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvp);
void dmkvvl(int nbasis, int kb, double *tn, int lang, struct bsp_basis *bss, struct vxaimp *pspa, double **pV, double **Vp, double **pVp);

double vpotloc(int z, double r, struct gauss_basis_set *gbs, struct orbital *orb);
double vpolhe(double r);
double param_polne(double r);
/*---------------------------------------------------------------------
   Evaluation of Bsplines overlap matrix ss by means of Gauss-Legendre
   integration.
---------------------------------------------------------------------*/

void mkss(int nbasis, int kb, double *tn, double **ss, struct bsp_basis *bss)
{
   double t,tmin,tmax,dt;
   double bibj,sum;
   int i,j,q;
   int ib,jb;

   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++){
      ib = bss[i].ib;					//associated bslines index for i
      jb = bss[j].ib;					//associated bslines index for j
      tmin = MAX(bss[i].tmin,bss[j].tmin);
      tmax = MIN(bss[i].tmax,bss[j].tmax);
      if(tmin >= tmax) ss[i][j] = 0.0;
      else{
         dt = tmax - tmin;
         sum = 0.0;
         for(q=0;q<QMAX_GAUSS;q++){
            t = (tmin + tmax + dt*egl[q])/2.0;
            bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
            sum += wgl[q]*bibj;}
         ss[i][j] = sum*dt/2.0;}
      ss[j][i] = ss[i][j];}

}



/*---------------------------------------------------------------------
   Evaluation of Bsplines kinetic energy matrix tt with additiional
   Bloch correction to ensures Hermicity. Gauss-Legendre integration
   is performed after analytitcal integration by part to remove thes
   outgoing flux.
   Implicit here is that the b/r term from Bloch operator assumes b = 0.
---------------------------------------------------------------------*/

void mkttbloch(int nbasis, int lang, int kb, double *tn, double **tt, struct bsp_basis *bss)
{
   double t,tmin,tmax,dt;
   double dbibj,bibj,sum;
   int i,j,q;
   int ib,jb;
   
   double ang2;

   ang2 = lang*(lang+1.0); 
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++){
      ib = bss[i].ib;					//associated bslines index for i
      jb = bss[j].ib;					//associated bslines index for j
      tmin = MAX(bss[i].tmin,bss[j].tmin);
      tmax = MIN(bss[i].tmax,bss[j].tmax);
      if(tmin >= tmax) tt[i][j] = 0.0;
      else{
         dt = tmax - tmin;
         sum = 0.0;
         for(q=0;q<QMAX_GAUSS;q++){
            t = (tmin + tmax + dt*egl[q])/2.0;
            bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
            dbibj = dbspline(tn,ib,kb,t)*dbspline(tn,jb,kb,t);
            sum += wgl[q]*dbibj;
            sum += wgl[q]*bibj*ang2/(t*t);}
         tt[i][j] = sum*dt/4.0;}			//extra factor 1/2 for kinetic energy
      tt[j][i] = tt[i][j];}
}






/*---------------------------------------------------------------------
   Evaluation of Bsplines potential matrix vv  for a local potential.
   ang2 = l*(l+1.0)/2.0;				//centrifugal potential
   Note that the centrifugal term is evaluated at the same time and
   added to the potential. Gauss-Legendre integration is performed.
---------------------------------------------------------------------*/
void mkvvl(int nbasis, int kb, double *tn, int lang, struct bsp_basis *bss, struct vxaimp *pspa, struct param_pola *ppola, double **pola, double **ionique, double **coulomb)
{
   char *atom;

   double t,tmin,tmax,dt;
   double ion, coul;
   double bibj, dbibj, ang2;
   double vp, vc, vi;
   double sum, sum2, sum3;
   int i,j,q;
   int ib,jb;
   int iprn = 0;
   FILE *fout;

   fout = fopen("pola.dat","w");
   sum = 0.0;
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
   {
      ib = bss[i].ib;					//associated bslines index for i
      jb = bss[j].ib;					//associated bslines index for j
      tmin = MAX(bss[i].tmin,bss[j].tmin);
      tmax = MIN(bss[i].tmax,bss[j].tmax);
      if(tmin >= tmax) pola[i][j] =  coulomb[i][j] =  ionique[i][j] = 0.0;
      else
      {
         dt = tmax - tmin;
         sum = sum2 = sum3 =  0.0;
         vc = vp = vi = 0.0;
         for(q=0;q<QMAX_GAUSS;q++)
         {
            t = (tmin + tmax + dt*egl[q])/2.0;
            bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
            vcpot(pspa,t,&ion,&coul);

            vp = polapot(ppola,t);
            if(POLA==1) fprintf(fout,"%lf \t %lf \n",t,vp);
            vi = ion;
            vc = coul;
            sum += wgl[q]*bibj*vp;
            sum2 += wgl[q]*bibj*vc;
            sum3 += wgl[q]*bibj*vi;
         }
         pola[i][j] = sum*dt/2.0;
         coulomb[i][j] = sum2*dt/2.0;
         ionique[i][j] = sum3*dt/2.0;
      }
      pola[j][i] = pola[i][j];
      coulomb[j][i] = coulomb[i][j];
      ionique[j][i] = ionique[i][j];
   }

   fclose(fout);

}


void dmkvvl(int nbasis, int kb, double *tn, int lang, struct bsp_basis *bss, struct vxaimp *pspa, double **pV, double **Vp, double **pVp)
{
   double t,tmin,tmax,dt;
   double ion, coul;
   double bi, bj, dbi, dbj;
   double vi,vc;
   double sum, sum2, sum3;
   int i,j,q;
   int ib,jb;
   int iprn = 0;

   sum = sum2 = sum3 = 0.0;

   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
   {
      ib = bss[i].ib;					//associated bslines index for i
      jb = bss[j].ib;					//associated bslines index for j
      tmin = MAX(bss[i].tmin,bss[j].tmin);
      tmax = MIN(bss[i].tmax,bss[j].tmax);
      if(tmin >= tmax) pV[i][j] = Vp[i][j] = pVp[i][j] = 0.0;
      else
      {
         dt = tmax - tmin;
         sum = sum2 = sum3 = 0.0;
         vi = vc = 0.0;
         for(q=0;q<QMAX_GAUSS;q++){
            t = (tmin + tmax + dt*egl[q])/2.0;
            bi = bspline(tn,jb,kb,t);
            bj = dbspline(tn,jb,kb,t);
            dbi = (dbspline(tn,ib,kb,t)*t - bspline(tn,ib,kb,t))/(t*t);
            dbj = (dbspline(tn,jb,kb,t)*t - bspline(tn,jb,kb,t))/(t*t);
            vcpot(pspa,t,&ion,&coul);

            vi = ion;
            vc = coul;
            sum += wgl[q]*(t*t*dbi*bj*(vi+vc));             // no minus since id/dr (p^+)
            sum2 += -wgl[q]*(t*t*bi*dbj*(vi+vc));           //minus come from -id/dr (p)
            sum3 += wgl[q]*(t*t*dbi*dbj*(vi+vc));}
         pV[i][j] = sum*dt/2.0;
         Vp[i][j] = sum2*dt/2.0;
         pVp[i][j] = sum3*dt/2.0;
      }
      pV[j][i] = pV[i][j];
      Vp[j][i] = Vp[i][j];
      pVp[j][i] = pVp[i][j];
   }
}

/*---------------------------------------------------------------------------
   bbvpa() compute the Pauli répulsion potential between 2 basis splines, for
   the whole basis set and returnsvoid hhaimp(int lang, int nb, int kb,  double *tn, double **ss, struct bsp_basis *bss, struct vxaimp *pspa, double **hh) the matrix vvp. The computation is done for
   a given angular momentum value l.
   The Pauli repulsion operator is simply the density matrix operator:
	vvp = Sum_l,ab |Gla> Pl,ab <Glb|
   where l is limited to the angular momenta of the occupied orbitals.
---------------------------------------------------------------------------*/
void mkvvpa(int nb, int kb,  double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvp)
{
   int lamax;
   struct vshell *shla;						//A copy of basis shell for occupied orbitals
   double **rho;						//A copy of density matrix
   double **wba;
   double **biga;
   double epauli;
//--------------
   int ia,ja,na;
   int ib,jb;
   double a;

   lamax = pspa->lbmax;
   for(ib=1;ib<=nb;ib++) for(jb=1;jb<=nb;jb++) vvp[ib][jb] = 0.0;
   printf("Pauli:   \t lamax = %d \t l = %d\n",lamax,l);

   if(l <= lamax)						//check whether l matches one of the occupied angular momentum
      {
      epauli = pspa->epauli;
      na = pspa->nlb[l];					//for bookkeeping
      shla = pspa->shlb[l];
      rho = pspa->rho[l];
      biga = dmatrix(1,nb,1,na);
      wba = dmatrix(1,nb,1,na);
//-------------------------------------
      for(ib=1;ib<=nb;ib++) for(ia=1;ia<=na;ia++){
         a = shla[ia].a;
         biga[ib][ia] = bspgto(kb,tn,bss[ib].ib,a,l);}		//bspline-RGTO overlappola
      for(jb=1;jb<=nb;jb++) for(ia=1;ia<=na;ia++){		//cumulate first half product
         wba[jb][ia] = 0.0;
         for(ja=1;ja<=na;ja++) wba[jb][ia] += biga[jb][ja]*rho[ia][ja];}
      for(ib=1;ib<=nb;ib++) for(jb=1;jb<=nb;jb++){		//cumulate second half product
         for(ia=1;ia<=na;ia++) vvp[ib][jb] += biga[ib][ia]*wba[jb][ia];
	 vvp[ib][jb] *= epauli;}
//-------------------------------------
      free_dmatrix(biga,1,nb,1,na);
      free_dmatrix(wba,1,nb,1,na);}
}

/*---------------------------------------------------------------------------
   bbvex() compute the exchange potential between 2 basis splines, for the
   whole basis set and returns the exchange matrix vvx. The computation is
   done for a given angular momentum value l.
   The exchange operator is assumed to be on the form:
	vvx = Sum_ab |Gla> Kl,ab <Glb|
---------------------------------------------------------------------------*/
void mkvvex(int nb, int kb, double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx)
{
   int lxmax;
   struct vshell *shlx;						//A copy of the extraction shell for the chosen "l"
   double **vvlx;						//A copy of the extracted exchange potential for "l"
   double **wbx;
   int ix,jx,nx;
   int ib,jb;
   double a;
   double **biga;

   lxmax = pspa->lxmax;
   for(ib=1;ib<=nb;ib++) for(jb=1;jb<=nb;jb++) vvx[ib][jb] = 0.0;
   printf("Exchange:\t lxmax = %d \t l = %d\n",lxmax,l);

   if(l <= lxmax)						//check whether l matches one of the occupied angular momentum
      {
      nx = pspa->nlx[l];
      shlx = pspa->shlx[l];
      vvlx = pspa->vvx[l];
      wbx = dmatrix(1,nb,1,nx);
//-------------------------------------
      biga = dmatrix(1,nb,1,nx);
      for(ib=1;ib<=nb;ib++) for(ix=1;ix<=nx;ix++)
         {
         a = shlx[ix].a;
         biga[ib][ix] = bspgto(kb,tn,bss[ib].ib,a,l);		//bspline-RGTO overlap
         }
      for(jb=1;jb<=nb;jb++) for(ix=1;ix<=nx;ix++)		//cumulate first half product
         {
         wbx[jb][ix] = 0.0;
         for(jx=1;jx<=nx;jx++) wbx[jb][ix] += biga[jb][jx]*vvlx[ix][jx];
         }
      for(ib=1;ib<=nb;ib++) for(jb=ib;jb<=nb;jb++)		//cumulate second half product
         {
         vvx[ib][jb] = 0.0;
         for(ix=1;ix<=nx;ix++) vvx[ib][jb] -= biga[ib][ix]*wbx[jb][ix];
         vvx[jb][ib] = vvx[ib][jb];				//set transposed element
         }
   //prmatxr2(11,20,vvx);
//-------------------------------------
      free_dmatrix(biga,1,nb,1,nx);
      free_dmatrix(wbx,1,nb,1,nx);
      }
}
