void build_potential(int zat, int nbasis, int nspline, int kb, int lang, double *tn, double **ss, double **tt, double **pola, double **ionique, double **coulomb, double **vvx, double **vvr,  double **vvc, struct bsp_basis *bss, struct bsp *bspline, struct basis_extraction *be, struct vxaimp *pspa, struct potcor *vcc, struct param_pola *ppola, double **dkh_hdk2, double **vcrelat);
void build_eigen(int l, int nbasis, double *eorb, double **corb, double **ss, double **tt, double **pola, double **ionique, double **coulomb, double **vvx, double **vvr, double **vvc, double **hh1, double **dkh_hdk2, double **vcrelat);

void build_potential(int zat, int nbasis, int nspline, int kb, int lang, double *tn, double **ss, double **tt, double **pola, double **ionique, double **coulomb, double **vvx, double **vvr,  double **vvc, struct bsp_basis *bss, struct bsp *bspline, struct basis_extraction *be, struct vxaimp *pspa, struct potcor *vcc, struct param_pola *ppola, double **dkh_hdk2, double **vcrelat)
{
    int i,j;
    static int init=0;

    if(init==0)
    {
        mkss(nbasis,kb,tn,ss,bss);
        init=1;
    }

    mkttbloch(nbasis,lang,kb,tn,tt,bss);
    mkvvl(nbasis,kb,tn,lang,bss,pspa,ppola,pola,ionique,coulomb);
    bbvexact(nbasis,kb,tn,lang,bss,pspa,vvx);
    mkvvpa(nbasis,kb,tn,lang,bss,pspa,vvr);

    mkvvcor(nbasis,nspline,kb,lang,tn,ss,vvc,bss,bspline,vcc,ppola); 

    DKH_construct(nbasis,kb,tn,lang,ss,tt,ionique,coulomb,dkh_hdk2,bss,pspa);
    //vdkh1(lang,zat,nbasis,kb,tn,vcrelat,bss);
}


void build_eigen(int l, int nbasis, double *eorb, double **corb, double **ss, double **tt,  double **pola, double **ionique, double **coulomb, double **vvx, double **vvr, double **vcc, double **hh1, double **dkh_hdk2, double **vcrelat)
{   
    int i,j;
    double **rel_corr;

    rel_corr = dmatrix(1,nbasis,1,nbasis);
    switch(DKH)
    {
        case 0 :
            for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++){hh1[i][j] = tt[i][j] + ionique[i][j] + coulomb[i][j] + pola[i][j] + vvx[i][j] + vvr[i][j] + vcc[i][j];}
        break;

        case 1 : 
            for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) rel_corr[i][j] = dkh_hdk2[i][j] - (tt[i][j] + ionique[i][j] + coulomb[i][j]);
            for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++){hh1[i][j] = tt[i][j] + ionique[i][j] + coulomb[i][j] + pola[i][j] + vvx[i][j] + vvr[i][j] + vcc[i][j] + rel_corr[i][j];}
        break;
        
        case 2 : 
            for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++){hh1[i][j] = tt[i][j] + ionique[i][j] + coulomb[i][j] + pola[i][j] + vvx[i][j] + vvr[i][j] + 0.0*vcc[i][j] + vcrelat[i][j];}
        break;
    }

    //for(i=1;i<=3;i++) printf("%lf\n",ss[i][i]);
    //getchar();
    diagoz(nbasis,ss,hh1,eorb,corb);

    printf("\nL = %d\n\n", l);
    for(i=1;i<=10;i++){printf("i=%d \t ei=%e \t de=%e \n",i,eorb[i],eorb[i]-eorb[1]);}
    printf("\n\n");
}
