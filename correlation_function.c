void init_vab(double **ss, struct potcor *vc, struct param_pola *ppola);
void mkvvcor(int nbasis, int nspline, int kb, int lang, double *tn, double **ss, double **vcor, struct bsp_basis *bss, struct bsp *bspl, struct potcor *vc, struct param_pola *ppola);

void init_vab(double **ss, struct potcor *vc, struct param_pola *ppola)
{
    char kw0[2];
    char str[2048];

    int i,j,k,l;
    int temp, dtemp;
    double ttemp;

    char *sf, *fname;
    FILE *fin;

    strcpy(kw0,"#");
    fname = (char *) malloc(200*sizeof(char));
    sprintf(fname,"./potcor/%s/potcor.in",ppola->atom);
    printf("correction file : %s\n", fname);

    if((fin = fopen(fname,"r")) == NULL)
    {
      printf("Cannot open input file : %s\n",fname);
      exit(1);
    }

    sf = fgets(str,511,fin);
    while(strncmp(str,kw0,1) == 0) sf = fgets(str,511,fin);
    sscanf(str,"%d",&temp);
    vc->lmax = temp;


    for(k=0;k<=vc->lmax;k++)
    {
        /*  Recuperation du LMAX + Taille de la Base    */
        sf = fgets(str,511,fin);
        while(strncmp(str,kw0,1) == 0) sf = fgets(str,511,fin);
        sscanf(str,"%d %d",&temp,&dtemp);
        l = temp;
        vc[l].size_be = dtemp;

        /*  Recuperation des arguments de la base    */

        for(i=1;i<=vc[l].size_be;i++)
        {
           fscanf(fin,"%lf   ",&ttemp);
           vc[l].arg[i] = ttemp;
        }

        /*  Recuperation du potentiel de correction fitté    */
        for(i=1;i<=vc[l].size_be;i++)
        {
            for(j=1;j<=vc[l].size_be;j++)
            {
                fscanf(fin,"%lf  ",&ttemp);
                vc[l].vab[i][j] = ttemp;
            }
        }            
    }

}

void mkvvcor(int nbasis, int nspline, int kb, int lang, double *tn, double **ss, double **vcor, struct bsp_basis *bss, struct bsp *bspl, struct potcor *vc, struct param_pola *ppola)
{
    static int initia=0;

    int q;
    int i,j,k,l;
    char *fname;
    
    double sum=0.0;
    double t,tmin,tmax,dt;
    double **temp, **ovlp;
 
    /* PRE-CALCUL DES VALEURS DES BSP POUR UN GAUSS-LEGENDRE */ 
    GL_Bsp(nspline,kb,tn,bspl);


    /* INITIALISATION DES VARIABLES DU POTENTIEL DE CORRECTION */
    
    init_vab(ss,vc,ppola);

    temp=dmatrix(1,nspline,1,nspline);
    ovlp=dmatrix(1,nspline,1,vc[lang].size_be);

    /* PRE-CALCUL DU RECOUVREMENT GAUSSIENNE-BSP POUR UN GAUSS-LEGENDRE */
    for(j=1;j<=vc[lang].size_be;j++){
        for(i=1;i<=nspline;i++){
            tmin = tn[i];
            tmax = tn[i+1+kb];
            dt = tmax - tmin;
            sum=0.0;
            for(q=0;q<QMAX_GAUSS;q++){
                t = (tmin + tmax + dt*egl[q])/2.0;
                sum += powint(t,lang)*bspl[i].BspValue[q]*exp(-vc[lang].arg[j]*t*t);}
            ovlp[i][j]=0.5*sum*dt;}}

    /*FACTORISATION DU TERME ANGULAIRE DU PRECEDENT CALCUL - GAIN CPU THEORIQUE */
    switch(lang){
        case 0 :
            for(j=1;j<=vc[lang].size_be;j++) for(i=1;i<=nspline;i++) ovlp[i][j] *= gtonorm(vc[lang].arg[j],0,0,0)*sqrt(4.*pi);
            break;
        
        case 1 :
            for(j=1;j<=vc[lang].size_be;j++) for(i=1;i<=nspline;i++) ovlp[i][j] *= gtonorm(vc[lang].arg[j],1,0,0)*sqrt(4.*pi/3.);
            break;

        case 2 :
            for(j=1;j<=vc[lang].size_be;j++) for(i=1;i<=nspline;i++) ovlp[i][j] *= gtonorm(vc[lang].arg[j],2,0,0)*sqrt(16.*pi/5.);
            break; 

        case 3 :
            for(j=1;j<=vc[lang].size_be;j++) for(i=1;i<=nspline;i++) ovlp[i][j] *= gtonorm(vc[lang].arg[j],3,0,0)*sqrt(16.*pi/7.);
            break;
            
        case 4 :
            for(j=1;j<=vc[lang].size_be;j++) for(i=1;i<=nspline;i++) ovlp[i][j] *= gtonorm(vc[lang].arg[j],3,0,0)*sqrt(16.*pi/7.);
            break;}    

    for(i=1;i<=nspline;i++){
        for(j=i;j<=nspline;j++){
            sum = 0.0;
            for(k=1;k<=vc[lang].size_be;k++){
                for(l=1;l<=vc[lang].size_be;l++){
                    sum += ovlp[i][k]*vc[lang].vab[k][l]*ovlp[j][l];}}
                    temp[i][j] = sum;
                    temp[j][i] = temp[i][j];}}

    for(i=1;i<=nbasis;i++){
        for(j=i;j<=nbasis;j++){
            vcor[i][j] = temp[bss[i].ib][bss[j].ib];
            vcor[j][i] = vcor[i][j];}}
}