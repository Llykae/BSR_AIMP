void build_PS(int l, int nbasis, double *E, double *eorb, double **corb, struct phase_shift *d0, struct bsp_basis *bss, struct vxaimp *pspa, struct param_pola *ppola);
void build_TCS(double *E, struct phase_shift *d0);
void build_DCS(struct param_pola *ppola, struct phase_shift *d0);
void build_DCS_bis(struct param_pola *ppola);

void build_PS(int l, int nbasis, double *E, double *eorb, double **corb, struct phase_shift *d0, struct bsp_basis *bss, struct vxaimp *pspa, struct param_pola *ppola)
{
    int i,n,ke,j;
    double zk; 
    double jl,nl,djl,dnl;
    double sl, dsl, cl, dcl;
    double aradius, rmtxe;

    double k=0.0;
    double e, dephs;
    double Etemp[Nk];

    int size;
    char trash[255], str[255];
    double *ee, *Ee;
    double Eigen;
    FILE *input;


    switch(TDSWC)
    {
      case 0 :     
        for(ke=1;ke<=Nk;ke++)
        {
            e = (1.0*ke*ke)/(Nk*Nk);
            E[ke] = e;
        
            aradius = bss[nbasis].tmax;			        //outer boundary
            bslsph(l,aradius*sqrt(2.0*e),&jl,&nl,&djl,&dnl);         //Get Neumann and Spherical Bessel func
            rmtxe = 0.0;
            for(n=1;n<=nbasis;n++){rmtxe += corb[nbasis][n]*corb[nbasis][n]/(eorb[n] - e);}
            
            rmtxe /= 2.0*aradius;
            zk = aradius*sqrt(2.0*e);
            

            sl = zk * jl;
            dsl = jl + zk*djl;
            cl = - zk * nl;
            dcl = - nl - zk*dnl;
            
            dephs = atan2(rmtxe*zk*dsl - sl,cl - rmtxe*zk*dcl);         //phase shift from inner space by R-matrix
            d0[l].ps[ke] = dcalogero(l,sqrt(2.0*e),dephs,aradius,pspa,ppola); //phase shift from outer space by Calogero

            if(d0[l].ps[ke] > pi/2.0) d0[l].ps[ke] -= pi;
            if(d0[l].ps[ke] <= -pi/2.0) d0[l].ps[ke] += pi;
        }
      break;

      case 1 :

          /* INITIALISATION USING INPUT : ./input/DCS_input.in */
          input=fopen("./input/DCS_input.in","r");
          fgets(trash,255,input);
          fgets(str,255,input);
          sscanf(str,"%d",&size);

          Ee = dvector(1,size);
          ee = dvector(1,size);

          fgets(trash,255,input);

          for(i=1;i<=size;i++){
            fgets(str,255,input);
            sscanf(str,"%lf  %lf",&(Ee[i]),&(ee[i]));}
          fclose(input);
          
          for(ke=1;ke<=size;ke++){

            Eigen = 0.5*Ee[ke]*Ee[ke];

            aradius = bss[nbasis].tmax;			        //outer boundary
            bslsph(l,aradius*Ee[ke],&jl,&nl,&djl,&dnl);         //Get Neumann and Spherical Bessel func
            rmtxe = 0.0;
            for(n=1;n<=nbasis;n++){rmtxe += corb[nbasis][n]*corb[nbasis][n]/(eorb[n] - Eigen);}
          
            rmtxe /= 2.0*aradius;
            zk = aradius*Ee[ke];
          

            sl = zk * jl;
            dsl = jl + zk*djl;
            cl = - zk * nl;
            dcl = - nl - zk*dnl;
            
            dephs = atan2(rmtxe*zk*dsl - sl,cl - rmtxe*zk*dcl);         //phase shift from inner space by R-matrix
            d0[l].ps[ke] = dcalogero(l,sqrt(2.0*e),dephs,aradius,pspa,ppola); //phase shift from outer space by Calogero
          
            if(d0[l].ps[ke] > pi/2.0) d0[l].ps[ke] -= pi;
            if(d0[l].ps[ke] <= -pi/2.0) d0[l].ps[ke] += pi;}
          break;
        }

}


void build_TCS(double *E, struct phase_shift *d0)
{
    int ke,l;           
    double sd, sum=0.0;
    double *sigma;

    FILE *fd;
    fd=fopen("./ScatteringResult/Scatteringtools.out","w");
  
    sigma=dvector(1,Nk);

    for(ke=1;ke<Nk;ke++){
      sum = 0.0;
      for(l=0;l<=LMAX;l++){
	      sd=sin(d0[l].ps[ke]);
	      sum+=(2.*l+1)*sd*sd;}

      sigma[ke]=2.0*pi/E[ke]*sum;                                         //total cross section
      printf("k=%d \t E=%e \t sigma=%e\n",ke,E[ke],sigma[ke]);

      fprintf(fd,"\n%lf",E[ke]);
      for(l=0;l<=LMAX;l++) fprintf(fd,"\t%e",d0[l].ps[ke]);
      fprintf(fd,"\t%e",sigma[ke]);}
    fclose(fd);
    
    free_dvector(sigma,1,Nk);
    
}

void build_DCS(struct param_pola *ppola, struct phase_shift *d0)
{
    int i,j,k;
    int size;

    double *E,*e;
    double sd,ssd,DCS;

    char *fname;
    char *eigen;
    char trash[255], str[255];

    FILE *input;
    FILE *output;

    /* INITIALISATION USING INPUT : ./input/DCS_input.in */
    input=fopen("./input/DCS_input.in","r");
    fgets(trash,255,input);
    fgets(str,255,input);
    sscanf(str,"%d",&size);

    E = dvector(1,size);
    e = dvector(1,size);

    fgets(trash,255,input);
    for(i=1;i<=size;i++){
      fgets(str,255,input);
      sscanf(str,"%lf  %lf",&(E[i]),&(e[i]));}
    fclose(input);


    /* CALCULATION OF THE DCS FOR EACH ENERGY OF THE INPUT*/
    fname = (char *) malloc(200*sizeof(char));
    eigen = (char *) malloc(100*sizeof(char));

    for(i=1;i<=size;i++){
      double_to_char(e[i],eigen);
      sprintf(fname,"./output/DCS/%s/DCS_%seV.out",ppola->atom,eigen);
      printf("%s\n",fname);
      getchar();
      output = fopen(fname,"w");

      for(j=1;j<=180;j++)
      {
        DCS = 0.0;
        sd = 0.0;
        ssd = 0.0;
        for(k=0;k<=LMAX;k++)
        {
          sd += (2.0*k+1.0)*sin(2.0*d0[k].ps[i])*plgndr(k,cos(pi*j/180.));
          ssd += (2.0*k+1.0)*(cos(2.0*d0[k].ps[i])-1.0)*plgndr(k,cos(pi*j/180.));
        }
        DCS = (0.25/(E[i]*E[i]))*(sd*sd+ssd*ssd);
        fprintf(output,"%lf   \t    %lf\n",1.0*j,DCS);
      }
      fclose(output);
    }  
}

void build_DCS_bis(struct param_pola *ppola)
{
    int i,n,k;
    int j = 0;
    int size=SIZE;
    int lmax=DCS_LMAX;

    double temp, temp2;
    double *E,*e;
    double sd,ssd,DCS;

    double **d0;
    char *fname, *finput, *fname_directory, *fname_input;
    char *eigen;
    char trash[255], str[255];

    FILE *input;
    FILE *output;

    finput = (char *) malloc(200*sizeof(char));;
    finput = EXTRACTION_DCS_FILE;
    
    /* INITIALISATION USING INPUT : ./input/DCS_input.in */
    fname_input = (char *) malloc(200*sizeof(char));
    sprintf(fname_input,"./data/electron/He/DCS/%s",finput);
    printf("%s\n",fname_input);
    input=fopen(fname_input,"r");

    d0 = dmatrix(1,size,1,lmax);
    E = dvector(1,size);
    e = dvector(1,size);

    for(i=1;i<=size;i++)
    {
      fscanf(input,"%lf \t %lf \t %lf \t %lf \t %lf\n", &(temp2), &(d0[i][1]), &(d0[i][2]), &(d0[i][3]), &temp);  
      E[i] = sqrt(temp2/13.605);
    }
    fclose(input);


    /* CALCULATION OF THE DCS FOR EACH ENERGY OF THE INPUT*/
    fname = (char *) malloc(200*sizeof(char));
    fname_directory = (char *) malloc(200*sizeof(char));
    eigen = (char *) malloc(100*sizeof(char));

    fname_directory = EXTRACTION_DCS_NAME;
    for(i=1;i<=size;i++){
      sprintf(fname,"./output/DCS/%s/%s/DCS_%lfeV.out",ppola->atom,fname_directory,(E[i]*E[i]*13.605));
      printf("%s\n",fname);
      getchar();
      output = fopen(fname,"w");

      j = 0;
      for(n=1;n<=45;n++)
      {
        j += 4;
        DCS = 0.0;
        sd = 0.0;
        ssd = 0.0;
        for(k=0;k<=LMAX;k++)
        {
          sd += (2.0*k+1.0)*sin(2.0*d0[i][k+1])*plgndr(k,cos(pi*j/180.));
          ssd += (2.0*k+1.0)*(cos(2.0*d0[i][k+1])-1.0)*plgndr(k,cos(pi*j/180.));
        }
        DCS = (0.25/(E[i]*E[i]))*(sd*sd+ssd*ssd);
        fprintf(output,"%lf   \t    %lf\n",1.0*j,DCS);
      }
      fclose(output);
    }  
}

void double_to_char(double f,char *buffer)
{
    gcvt(f,15,buffer);
}