void GL_Bsp(int nspline, int kb, double *tn, struct bsp *b);
void GL_GB(int swc, int k, int nspline, int kb, double *tn, double **ovrlp, struct bsp *b, struct vxaimp *pspa);

void GL_Bsp(int nspline, int kb, double *tn, struct bsp *b)
{
  int i,q;
  double tmin, tmax, t, dt;

  for(i=1;i<=nspline;i++){
    tmin = tn[i];
    tmax = tn[i+1+kb];
    dt = tmax - tmin;
    for(q=0;q<QMAX_GAUSS;q++){
      t = (tmin + tmax + dt*egl[q])/2.0;
      b[i].BspValue[q] = wgl[q]*t*bspline(tn,i,kb,t);}}
}

void GL_GB(int lang, int k, int nspline, int kb, double *tn, double **ovrlp, struct bsp *b, struct vxaimp *pspa)
{
    int i,q;
    double t,tmin,tmax,dt;
    double sum=0.0;
    struct vshell *shlx;

    switch(lang){

        case 0 :
            shlx=pspa->shlx[lang];
            for(i=1;i<=nspline;i++){
                tmin = tn[i];
                tmax = tn[i+1+kb];
                dt = tmax - tmin;
                sum=0.0;
                for(q=0;q<QMAX_GAUSS;q++){
                    t = (tmin + tmax + dt*egl[q])/2.0;
                    sum += b[i].BspValue[q]*exp(-shlx[k].a*t*t);}
                ovrlp[i][k]=0.5*gtonorm(shlx[k].a,0,0,0)*sqrt(4*pi)*sum*dt;}   
        break;

        case 1 :
            shlx=pspa->shlx[lang];
            for(i=1;i<=nspline;i++){
                tmin = tn[i];
                tmax = tn[i+1+kb];
                dt = tmax - tmin;
                sum=0.0;
                for(q=0;q<QMAX_GAUSS;q++){
                    t = (tmin + tmax + dt*egl[q])/2.0;
                    sum += t*b[i].BspValue[q]*exp(-shlx[k].a*t*t);}
                ovrlp[i][k]=0.5*gtonorm(shlx[k].a,1,0,0)*sqrt(4.*pi/3.)*sum*dt;}
        break;
        
        case 2 : 
            shlx=pspa->shlx[lang];
            for(i=1;i<=nspline;i++){
                tmin = tn[i];
                tmax = tn[i+1+kb];
                dt = tmax - tmin;
                sum=0.0;
                for(q=0;q<QMAX_GAUSS;q++){
                    t = (tmin + tmax + dt*egl[q])/2.0;
                    sum += t*t*b[i].BspValue[q]*exp(-shlx[k].a*t*t);}
                ovrlp[i][k]=0.5*gtonorm(shlx[k].a,2,0,0)*sqrt(16.*pi/5.)*sum*dt;}
        break;}
}