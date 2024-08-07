void init_basis_extraction(char *fname, struct basis_extraction *be);
void read_Vx(char *fname, struct basis_extraction *be);
void read_basis_extraction(struct basis_extraction *be);
void read_Vpauli(double **density);





void read_Vx(char *fname, struct basis_extraction *be)
{
    int i,j;
    int nLine = 0;
    int pbl;
    char trash[255];
    FILE *fin;

    fin = fopen(fname, "r");
    nLine = found_occu("tvv2ex", fin);
    rewind(fin);

    for(i=1;i<=nLine+1;i++){fgets(trash,255,fin);}

    for(i=1;i<=be[0].ns;i++){
        for(j=1;j<=be[0].ns;j++){fscanf(fin,"%lf ",&(be[0].vs[i][j]));}}


    fgets(trash,255,fin);
    for(i=1;i<=be[0].np;i++){
        for(j=1;j<=be[0].np;j++){fscanf(fin,"%lf ",&(be[0].vp[i][j]));}}


    fgets(trash,255,fin);
    for(i=1;i<=be[0].nd;i++){
        for(j=1;j<=be[0].nd;j++){fscanf(fin,"%lf ",&(be[0].vd[i][j]));}}

    fclose(fin);
}

void init_basis_extraction(char *fname, struct basis_extraction *be)
{
    int i,t;
    int nLine = 0;
    char trash[255];
    FILE *fin;

    fin = fopen(fname, "r");
    nLine = found_occu("tvv2ex", fin);
    rewind(fin);

    for(i=1;i<=nLine;i++){fgets(trash, 255, fin);}

    fscanf(fin,"%d %d",&t,&(be[0].ns));
    for(i=1;i<=be[0].ns+1;i++){fgets(trash,255,fin);}

    fscanf(fin,"%d %d",&t,&(be[0].np));
    for(i=1;i<=be[0].np+1;i++){fgets(trash,255,fin);}


    fscanf(fin,"%d %d",&t,&(be[0].nd));

    rewind(fin);
    nLine = found_occu("lambda", fin);
    rewind(fin);

    for(i=1;i<=nLine;i++){fgets(trash, 255, fin);};
    for(i=1;i<=be[0].ns;i++){fscanf(fin,"%d %lf",&t,&(be[0].args[i]));}
    for(i=1;i<=be[0].np;i++){fscanf(fin,"%d %lf",&t,&(be[0].argp[i]));}
    for(i=1;i<=be[0].nd;i++){fscanf(fin,"%d %lf",&t,&(be[0].argd[i]));}

    fclose(fin);
}

void read_Vpauli(double **density)
{
    int i,j;
    int nLine = 0;
    int pbl, lim;
    char trash[255], *fname;
    FILE *fin;

    fname = "vaimp.out_He6s6p2d1f";

    fin = fopen(fname, "r");
    nLine = found_occu("vpauli", fin);
    rewind(fin);

    for(i=1;i<=nLine;i++){fgets(trash,255,fin);}
    fscanf(fin,"%d %d", &pbl, &lim);

    for(i=1;i<=lim;i++){
        for(j=1;j<=lim;j++){fscanf(fin,"%lf ",&(density[i][j]));}}

    fclose(fin);
}

void read_basis_extraction(struct basis_extraction *be)
{
    char *fname;

    fname = "./vaimp.out_He6s6p2d1f";
    init_basis_extraction(fname,be);
    read_Vx(fname,be);
}