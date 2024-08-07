void read_input_pola(struct param_pola *ppola)
{
    int i,j;
    char filename[128];
    int itemp;
    double dtemp;
    FILE *fin;

    char kwin[15];
    char kwdipf1[15];
    char kwdipf2[15];
    char kwquaf1[15];
    char kwquaf2[15];
    char kwquaf3[15];
    char kwoctuf1[15];

    char *sf;
    char string[2048];

    strcpy(kwin,"# INPUT");
    strcpy(kwdipf1,"# DIPOLAR F1");
    strcpy(kwdipf2,"# DIPOLAR F2");
    strcpy(kwquaf1,"# QUADRI F1");
    strcpy(kwquaf2,"# QUADRI F2");
    strcpy(kwquaf3,"# QUADRI F3");
    strcpy(kwoctuf1,"# OCTU F1");

    sprintf(filename,"./input_pola/pola_%s.dat",ppola->atom);
    printf("pola file : %s \n",filename);

    /*  READ ATOMIC PARAMETER */
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwin,7) != 0) sf = fgets(string,511,fin);
    for(i=1;i<=4;i++)
    {
        fscanf(fin,"%lf\n",&dtemp);
        ppola->aab[i] = dtemp;
    }
    //printf("alphad : %+le \t alphaq : %+le \t beta : %+le\n\n",ppola->aab[1],ppola->aab[2],ppola->aab[3]);
    fclose(fin);


    /* READ DIPOLAR F1 */
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwdipf1,12) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->ndipol[1] = itemp;
    
    ppola->dipow[1] = ivector(1,ppola->ndipol[1]);
    ppola->dipg[1] = dvector(1,ppola->ndipol[1]);
    ppola->dipc[1] = dvector(1,ppola->ndipol[1]);
    for(i=1;i<=ppola->ndipol[1];i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->dipow[1][i]), &(ppola->dipg[1][i]), &(ppola->dipc[1][i]));
    //for(i=1;i<=ppola->ndipol[1];i++) printf("%d \t %lf \t %lf\n",ppola->dipow[1][i],ppola->dipg[1][i],ppola->dipc[1][i]);
    //printf("\n");
    fclose(fin);



    /*  READ DIPOLAR F2*/
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwdipf2,12) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->ndipol[2] = itemp;

    ppola->dipow[2] = ivector(1,ppola->ndipol[2]);
    ppola->dipg[2] = dvector(1,ppola->ndipol[2]);
    ppola->dipc[2] = dvector(1,ppola->ndipol[2]);
    for(i=1;i<=ppola->ndipol[2];i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->dipow[2][i]), &(ppola->dipg[2][i]), &(ppola->dipc[2][i]));
    //for(i=1;i<=ppola->ndipol[2];i++) printf("%d \t %lf \t %lf\n",ppola->dipow[2][i],ppola->dipg[2][i],ppola->dipc[2][i]);
    //printf("\n");
    fclose(fin);



    /*  READ QUADRI F1*/
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwquaf1,11) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->nquad[1] = itemp;

    ppola->quadpow[1] = ivector(1,ppola->nquad[1]);
    ppola->quadg[1] = dvector(1,ppola->nquad[1]);
    ppola->quadc[1] = dvector(1,ppola->nquad[1]);
    for(i=1;i<=ppola->nquad[1];i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->quadpow[1][i]), &(ppola->quadg[1][i]), &(ppola->quadc[1][i]));
    //for(i=1;i<=ppola->nquad[1];i++) printf("%d \t %lf \t %lf\n",ppola->quadpow[1][i],ppola->quadg[1][i],ppola->quadc[1][i]);
    //printf("\n");
    fclose(fin);




    /*  READ QUADRI F2*/
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwquaf2,11) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->nquad[2] = itemp;

    ppola->quadpow[2] = ivector(1,ppola->nquad[2]);
    ppola->quadg[2] = dvector(1,ppola->nquad[2]);
    ppola->quadc[2] = dvector(1,ppola->nquad[2]);
    for(i=1;i<=ppola->nquad[2];i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->quadpow[2][i]), &(ppola->quadg[2][i]), &(ppola->quadc[2][i]));
    //for(i=1;i<=ppola->nquad[2];i++) printf("%d \t %lf \t %lf\n",ppola->quadpow[2][i],ppola->quadg[2][i],ppola->quadc[2][i]);
    //printf("\n");
    fclose(fin);




    /*  READ QUADRI F3*/
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwquaf3,11) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->nquad[3] = itemp;

    ppola->quadpow[3] = ivector(1,ppola->nquad[3]);
    ppola->quadg[3] = dvector(1,ppola->nquad[3]);
    ppola->quadc[3] = dvector(1,ppola->nquad[3]);
    for(i=1;i<=ppola->nquad[3];i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->quadpow[3][i]), &(ppola->quadg[3][i]), &(ppola->quadc[3][i]));
    //for(i=1;i<=ppola->nquad[3];i++) printf("%d \t %lf \t %lf\n",ppola->quadpow[3][i],ppola->quadg[3][i],ppola->quadc[3][i]);
    //printf("\n");
    fclose(fin);


    /*  READ OCTU F1*/
    fin = fopen(filename, "r");
    if((fin = fopen(filename,"r")) == NULL)
    {
        printf("Cannot open input file : pola.dat\n");
        exit(1);
    }

    sf = fgets(string,511,fin);
    while(strncmp(string,kwoctuf1,9) != 0) sf = fgets(string,511,fin);
    fscanf(fin,"%d\n",&itemp);
    ppola->noctu = itemp;
    ppola->octupow = ivector(1,ppola->noctu);
    ppola->octug = dvector(1,ppola->noctu);
    ppola->octuc = dvector(1,ppola->noctu);
    for(i=1;i<=ppola->noctu;i++) fscanf(fin,"%d \t %lf \t %lf\n", &(ppola->octupow[i]), &(ppola->octug[i]), &(ppola->octuc[i]));
    fclose(fin);

}





double dipole(struct param_pola *ppola, double r)
{
    int i;
    double f1, f2;

    f1 = 1.0;
    for(i=1;i<=ppola->ndipol[1];i++) f1 -= pow(r,ppola->dipow[1][i])*ppola->dipc[1][i]*exp(-ppola->dipg[1][i]*r*r);
    
    f2 = 0.0;
    for(i=1;i<=ppola->ndipol[2];i++) f2 += pow(r,ppola->dipow[2][i])*ppola->dipc[2][i]*exp(-ppola->dipg[2][i]*r*r);

    return -0.5*ppola->aab[1]*(f1*f1/(r*r*r*r) + f2);
}

double quadripole(struct param_pola *ppola, double r)
{
    int i;
    double f1, f2, f3;

    f1 = 1.0;
    for(i=1;i<=ppola->nquad[1];i++) f1 -= pow(r,ppola->quadpow[1][i])*ppola->quadc[1][i]*exp(-ppola->quadg[1][i]*r*r);

    f2 = 0.0;
    for(i=1;i<=ppola->nquad[2];i++) f2 += pow(r,ppola->quadpow[2][i])*ppola->quadc[2][i]*exp(-ppola->quadg[2][i]*r*r);

    f3 = 0.0;
    for(i=1;i<=ppola->nquad[3];i++) f3 += pow(r,ppola->quadpow[3][i])*ppola->quadc[3][i]*exp(-ppola->quadg[3][i]*r*r);
    
    return -0.5*(ppola->aab[2]-6.*ppola->aab[3])*(f1*f1/(r*r*r*r*r*r)) - 0.5*ppola->aab[1]*f3 + 3.0*ppola->aab[3]*f2;
}

double octupole(struct param_pola *ppola, double r)
{
    int i;
    double f1;

    f1 = 0.0;
    for(i=1;i<=ppola->noctu;i++) f1 += pow(r,ppola->octupow[i])*ppola->octuc[i]*exp(-ppola->octug[i]*r*r);

    return 3.*ppola->aab[4]*f1;
}

double polapot(struct param_pola *ppola, double r)
{
    return dipole(ppola,r) + quadripole(ppola,r) + 0.0*octupole(ppola,r);
}