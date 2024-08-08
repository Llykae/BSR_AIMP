void readvaipp(struct vxaimp *pspa,struct param_pola *ppola);
//---
void bbvex(int nb, int kb,  double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx);
//---
void bbvcl(int nb, int kb,  double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvx);
double bspvcl(int kb,  double *tn, int ib, int jb, struct vxaimp *pspa, int l);
//double vcpot(struct vxaimp *pspa, double r);
double uclfunc(int l, double q, double r2);
//---
void bbvpa(int nb, int kb,  double *tn, int l, struct bsp_basis *bss, struct vxaimp *pspa, double **vvp);
//void hhaimp(int lang, int nb, int kb,  double *tn, double **ss, struct bsp_basis *bss, struct vxaimp *pspa, double **hh);
/*---------------------------------------------------------------------------
   readvaipp() reads the parameters of the ab initio model potential from
   the input file "genaimp/output/vaimp.out".
   The parameters are stored in the structure pspa.
---------------------------------------------------------------------------*/
void readvaipp(struct vxaimp *pspa,struct param_pola *ppola)
{
   char kw0[2];
   char kwvex[10];
   char kwden[8];
   char kwato[5];
   char kwepa[7];
//--------------
   char *fname;
   FILE *fin;
   char string[512];
   //char word[32];
   char *sf;
   //int ls,lw;
//--------------
   int ib,jb;
   int lb,lbmax,nlb;
//--------------
   int ix,jx;
   int lx,lxmax,nlx;
   //double ax;
   int itmp;
   double dtmp;
//--------------
   int iprint = 1;

   strcpy(kw0,"#");						//# = comment line indicator
   strcpy(kwato,"ATOM");					//ATOM keyword
   strcpy(kwden,"DENSITY");					//DENSITY keyword
   strcpy(kwepa,"EPAULI");					//EPAULI keyword
   strcpy(kwvex,"VEXCHANGE");					//VEXCHANGE keyword

//-------------------------------------
   fname = (char *) malloc(200*sizeof(char));
   sprintf(fname,"genaimp/output/vaimp_%s.out",ppola->atom);
   if((fin = fopen(fname,"r")) == NULL)
      {
      printf("Cannot open input file : genaimp/output/vaimp.out\n");
      exit(1);
      }
   sf = fgets(string,511,fin);
   while(strncmp(string,kwato,4) != 0) sf = fgets(string,511,fin);
   printf("%s\n",string);
   fscanf(fin,"%le",&dtmp);
   pspa->zi = dtmp;
   fscanf(fin,"%le",&dtmp);
   pspa->qe = dtmp;
   pspa->zc = pspa->zi - pspa->qe;
   printf("Zi=%+1.3f \t Qe=%+1.3f \t Zc=%+1.3f\n\n",pspa->zi,-pspa->qe,pspa->zc);
   fclose(fin);

//-------------------------------------
   if((fin = fopen(fname,"r")) == NULL)
      {
      printf("Cannot open input file : genaimp/output/vaimp.out\n");
      exit(1);
      }
   sf = fgets(string,511,fin);
   while(strncmp(string,kwden,6) != 0) sf = fgets(string,511,fin);
   printf("%s\n",string);
   sf = fgets(string,511,fin);
   while(strncmp(string,kw0,1) == 0) sf = fgets(string,511,fin);
   sscanf(string,"%d",&lbmax);
   printf("lbmax = %d \n",lbmax);
   pspa->lbmax = lbmax;
   for(lb=0;lb<=lbmax;lb++)
      {
      sf = fgets(string,511,fin);
      while(strncmp(string,kw0,1) == 0) sf = fgets(string,511,fin);
//printf("S: %s\n",string);
      sscanf(string,"%d %d",&itmp,&nlb);
      printf("lb=%d \t nlb=%d \n",lb,nlb);
      pspa->nlb[lb] = nlb;
      pspa->shlb[lb] = malloc((nlb+1)*sizeof(struct vshell));;
      pspa->rho[lb] = dmatrix(1,nlb,1,nlb);
      for(ib=1;ib<=nlb;ib++)
         {
         fscanf(fin,"%le",&dtmp);
         pspa->shlb[lb][ib].a = dtmp;
         pspa->shlb[lb][ib].l = lb;
         pspa->shlb[lb][ib].norm = rgtonorm(pspa->shlb[lb][ib].a,lb);
         if(iprint)  printf("lb=%d \t ab=%+1.4le \t normb=%+le \n",pspa->shlb[lb][ib].l,pspa->shlb[lb][ib].a,pspa->shlb[lb][ib].norm);
         }
      for(ib=1;ib<=nlb;ib++) for(jb=1;jb<=nlb;jb++)
         {
         fscanf(fin,"%le",&dtmp);
         pspa->rho[lb][ib][jb] = dtmp;
         }
      sf = fgets(string,511,fin);				//read end of last line including "\0"
      if(iprint) prmatxr2(1,nlb,pspa->rho[lb]);
//      getchar();
      rhonorm(lb,nlb,pspa->shlb[lb],pspa->rho[lb]);
      }
   fclose(fin);

//-------------------------------------
   if((fin = fopen(fname,"r")) == NULL)
      {
      printf("Cannot open input file : genaimp/output/vaimp.out\n");
      exit(1);
      }
   sf = fgets(string,511,fin);
   while(strncmp(string,kwepa,5) != 0) sf = fgets(string,511,fin);
   printf("%s\n",string);
   sf = fgets(string,511,fin);
   while(strncmp(string,kw0,1) == 0) sf = fgets(string,511,fin);
   sscanf(string,"%le",&dtmp);
   pspa->epauli = dtmp;
   printf("epauli = %+le \n",pspa->epauli);

//-------------------------------------
   if((fin = fopen(fname,"r")) == NULL)
      {
      printf("Cannot open input file : genaimp/output/vaimp.out\n");
      exit(1);
      }
   sf = fgets(string,511,fin);
   while(strncmp(string,kwvex,8) != 0) sf = fgets(string,511,fin);
   printf("%s\n",string);
   sf = fgets(string,511,fin);
   while(strncmp(string,kw0,1) == 0) sf = fgets(string,511,fin);
   sscanf(string,"%d",&lxmax);
   printf("lxmax = %d \n",lxmax);
   pspa->lxmax = lxmax;
   for(lx=0;lx<=lxmax;lx++)
      {
      sf = fgets(string,511,fin);
      while(strncmp(string,kw0,1) == 0) sf = fgets(string,511,fin);
//printf("S: %s\n",string);
      sscanf(string,"%d %d",&itmp,&nlx);
      printf("lx=%d \t nlx=%d \n",lx,nlx);
      pspa->nlx[lx] = nlx;
      pspa->shlx[lx] = malloc((nlx+1)*sizeof(struct vshell));;
      pspa->vvx[lx] = dmatrix(1,nlx,1,nlx);
      for(ix=1;ix<=nlx;ix++)
         {
         fscanf(fin,"%le",&dtmp);
         pspa->shlx[lx][ix].a = dtmp;
         pspa->shlx[lx][ix].l = lx;
         pspa->shlx[lx][ix].norm = rgtonorm(pspa->shlx[lx][ix].a,lx);
         if(iprint)  printf("lx=%d \t ax=%+1.4le \t normx=%+le \n",pspa->shlx[lx][ix].l,pspa->shlx[lx][ix].a,pspa->shlx[lx][ix].norm);
         }
      for(ix=1;ix<=nlx;ix++) for(jx=1;jx<=nlx;jx++)
         {
         fscanf(fin,"%le",&dtmp);
         pspa->vvx[lx][ix][jx] = dtmp;
         }
      sf = fgets(string,511,fin);				//read end of last line including "\0"
      if(iprint) prmatxr2(1,nlx,pspa->vvx[lx]);
//      getchar();
      }
   fclose(fin);

//--------------
}