
void sstransform(int nbasis, int kb, double **ss, double **mm, double **lmmr);
void aatransform(int nbasis, double **aa, double **mm, double **lmmr);

/*-----------------------------------------------------------------
   Relativistic correction for kinetic and ionic energy
   Cf. Nakajima et Hirao, Chemical Reviews 112, (2012) p. 385
   Classic:		Tc = p^2/2
   			Vc = -Z/r
   relativistic:	E0 = sqrt(p^2 c^2 + c^4)
   			Tr = E0 - c^2
			A = sqrt((E0 + c^2)/2E0)
			Vr = A.Vc.A
   Difference:		H = Tr + Vr - Tc - Vc
-----------------------------------------------------------------*/
void vrelc(int lang, double zat, int nbasis, int kb, double *tn, double **hh, struct bsp_basis *bss)
{
    double clight2;
    double **pp0,**ss,**vv;
    double **ttr,**vvr;
    double **cc0;
    double **aa,**kk,**ks,**cvvc,**avva;
    double *e0,*er,*e1,*em,*ae;
    double t,tmin,tmax,dt;
    double dbidbj,bibj,sum,sums,sumv;
    double lang2;
    int i,j,q;
    int n;
    int ib,jb;
    //---
    int iprint=1;

    ss = dmatrix(1,nbasis,1,nbasis);			//allocate overlap matrix
    pp0 = dmatrix(1,nbasis,1,nbasis);			//allocate one-electron P^2 matrix
    ttr = dmatrix(1,nbasis,1,nbasis);			//allocate relativistic kinetic energy matrix
    cc0 = dmatrix(1,nbasis,1,nbasis);			//allocate orbital coefficient matrix
    vv = dmatrix(1,nbasis,1,nbasis);			//allocate classical potential energy matrix
    vvr = dmatrix(1,nbasis,1,nbasis);			//allocate relativistic potential energy matrix
    avva = dmatrix(1,nbasis,1,nbasis);
    cvvc = dmatrix(1,nbasis,1,nbasis);
    aa = dmatrix(1,nbasis,1,nbasis);			//allocate intermediate matrix
    kk = dmatrix(1,nbasis,1,nbasis);
    ks = dmatrix(1,nbasis,1,nbasis);
    e0 = dvector(1,nbasis);
    em = dvector(1,nbasis);
    er = dvector(1,nbasis);
    ae = dvector(1,nbasis);

    clight2 = 137.036;
    clight2 *= clight2;
    lang2 = lang*(lang+1.0);
    for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
    {
        ib = bss[i].ib;                                   //associated bsplines index for i
        jb = bss[j].ib;                                   //associated bsplines index for j
        tmin = MAX(bss[i].tmin,bss[j].tmin);
        tmax = MIN(bss[i].tmax,bss[j].tmax);
        if(tmin >= tmax) pp0[i][j] = ss[i][j] = vv[i][j] = 0.0;
        else
        {
            dt = tmax - tmin;
            sum = sums = sumv = 0.0;
            for(q=0;q<QMAX_GAUSS;q++)
            {
                t = (tmin + tmax + dt*egl[q])/2.0;
                dbidbj = dbspline(tn,ib,kb,t)*dbspline(tn,jb,kb,t);
                bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
                sum += wgl[q]*(dbidbj + bibj*lang2/(t*t));
                sums += wgl[q]*bibj;
                sumv -= wgl[q]*bibj*zat/t;
            }
            pp0[i][j] = sum*dt/2.0;
            ss[i][j] = sums*dt/2.0;
            vv[i][j] = sumv*dt/2.0;
        }
        pp0[j][i] = pp0[i][j];					//classical momentum square in Bspline basis: P^2 (T0 = P^2/2
        ss[j][i] = ss[i][j];					//overlap matrix in Bspline basis
        vv[j][i] = vv[i][j];					//potential energy matrix in Bspline basis
    }


//-------------------------
    diagoz(nbasis,ss,pp0,e0,cc0);				//get P^2 basis: cc0
    
    //printf("lang = %d \t zat = %+1.3f \t kb=%d \n",lang,zat,kb);
    for(n=1;n<=nbasis;n++)
    {
        em[n] = sqrt(e0[n]*clight2 + clight2*clight2);		//relativistic free energy
        er[n] = em[n] - clight2;					//relativistic kinetic energy
        ae[n] = sqrt(0.5*(em[n]+clight2)/em[n]);			//A factor
        //if(n%10 == 0) printf("n=%d \t t0=%+le \t tr=%+le \n",n,e0[n]/2.0,er[n]);
    }

    //-------------------------
    for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)		//set relativistic kinetic energy in P^2 basis
    {
        kk[i][j] = aa[i][j] = 0.0;
    //      for(n=1;n<=nbasis;n++)  kk[i][j] += cc0[i][n]*cc0[j][n]*e0[n]/2.0;	//classical energy
        for(n=1;n<=nbasis;n++)
        {
            kk[i][j] += cc0[i][n]*cc0[j][n]*er[n];			//relativistic energy
            aa[i][j] += cc0[i][n]*cc0[j][n]*ae[n];			//ae in Bspline basis
        }
        kk[j][i] = kk[i][j];
        aa[j][i] = aa[i][j];
    //if(i<=20 && j<=(i+12)) printf("%d  %d\t %+le\n",i,j,ss[i][j]);
    }
    //printf("end K\n");
    sstransform(nbasis,kb,ss,kk,ttr);				//back in bsplines basis

    //-------------------------
   aatransform(nbasis,aa,vv,avva);				//A.V.A in Bspline basis
   sstransform(nbasis,kb,ss,avva,vvr);				//back in bsplines basis

//-------------------------
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
    {
      hh[i][j] = ttr[i][j] + vvr[i][j] - 0.5*pp0[i][j] - vv[i][j];	//difference Hamiltonina
      hh[j][i] = hh[i][j];
    }

    for(i=1;i<=10;i++) printf("%+le \t  %+le\n",0.5*pp0[i][i],ttr[i][i]);
//-------------------------
    /*if(iprint)
    {
        e1 = dvector(1,nbasis);
        diagoz(nbasis,ss,ttr,e1,cc0);
        for(n=1;n<=nbasis;n+=10)
        {
    //         printf("n=%d \t \t ae=%+le \t t0=%+le \t tr=%+le \n",n,ae[n],e0[n]/2.0,er[n]);
            printf("n=%d \t \t t0=%+le \t tr=%+le \t vv0=%+le \t vvr=%+le \t hh=%+le\n",n,0.5*pp0[n][n],ttr[n][n],vv[n][n],vvr[n][n],hh[n][n]);
        }
        free_dvector(e1,1,nbasis);
    }*/

    //---
    free_dmatrix(ss,1,nbasis,1,nbasis);
    free_dmatrix(pp0,1,nbasis,1,nbasis);
    free_dmatrix(ttr,1,nbasis,1,nbasis);
    free_dmatrix(vv,1,nbasis,1,nbasis);
    free_dmatrix(vvr,1,nbasis,1,nbasis);
    free_dmatrix(avva,1,nbasis,1,nbasis);
    free_dmatrix(cvvc,1,nbasis,1,nbasis);
    free_dmatrix(cc0,1,nbasis,1,nbasis);
    free_dmatrix(aa,1,nbasis,1,nbasis);
    free_dmatrix(kk,1,nbasis,1,nbasis);
    free_dmatrix(ks,1,nbasis,1,nbasis);
    free_dvector(e0,1,nbasis);
    free_dvector(em,1,nbasis);
    free_dvector(er,1,nbasis);
    free_dvector(ae,1,nbasis);
}


/*-----------------------------------------------------------------
   Relativistic correction for kinetic and ionic energy
   Cf. Nakajima et Hirao, Chemical Reviews 112, (2012) p. 385
   For the Darwin term, only the first basis function contributes by
   its derivative at r = 0. We use the defitinition given by Bethe
   and Salpeter page 59, formula 13.6, after symetrization of the
   operator:

   	Darwin -> Z/8c^2 . [db1/dr]^2(r=0)

   The comlpete term for DKH1 is:
   	DKH1 -> A.R.V.R.A - A.B.P.V.P.B.A

   Classic:		Tc = p^2/2
   			Vc = -Z/r
   relativistic:	E0 = sqrt(p^2 c^2 + c^4)
   			Tr = E0 - c^2
			A = sqrt((E0 + c^2)/2E0)
			B = c/(E0 + c^2)
			Vr = A.(Vc + R.Vc.R).A
   Difference:		H = Tr + Vr - Tc - Vc

   For moderate kinetic energy, A is of the order of unity and
   B is of the order of 1/c. The DKH1 expression is correct to all
   orders unlike the Pauli approximation, which assumes A = 1 and
   B = 1/2, and keeps only E0 = p^4/c^2 and the Darwin term limited
   to its expression for S waves. For L > 0, the term R.Vc.R brings
   a small corection.
-----------------------------------------------------------------*/
void vdkh1(int lang, double zat, int nbasis, int kb, double *tn, double **hh, struct bsp_basis *bss)
{
   double clight,clight2;
   double **pp0,**ss,**vv,**pvvp;
   double **ttr,**vvr,**vvd;
   double **cc0;
   double **aa,**kk,**ks,**avva;
   double **yy,**zz;
   double *e0,*er,*e1,*em,*ae,*be,*ab1;
   double t,tmin,tmax,dt;
   double dbidbj,bibj;
   double sum,sums,sumv;
   double lang2;
   double darwin;
   int i,j,q;
   int m,n;
   int ib,jb;
//---
   int iprint=1;

   ss = dmatrix(1,nbasis,1,nbasis);			//allocate overlap matrix
   pp0 = dmatrix(1,nbasis,1,nbasis);			//allocate one-electron P^2 matrix
   ttr = dmatrix(1,nbasis,1,nbasis);			//allocate relativistic kinetic energy matrix
   cc0 = dmatrix(1,nbasis,1,nbasis);			//allocate orbital coefficient matrix
   vv = dmatrix(1,nbasis,1,nbasis);			//allocate classical potential energy matrix
   pvvp = dmatrix(1,nbasis,1,nbasis);
   vvr = dmatrix(1,nbasis,1,nbasis);			//allocate relativistic potential energy matrix
   vvd = dmatrix(1,nbasis,1,nbasis);			//allocate relativistic Darwin potential energy matrix
   avva = dmatrix(1,nbasis,1,nbasis);
   aa = dmatrix(1,nbasis,1,nbasis);			//allocate intermediate matrix
   kk = dmatrix(1,nbasis,1,nbasis);
   ks = dmatrix(1,nbasis,1,nbasis);
   yy = dmatrix(1,nbasis,1,nbasis);
   zz = dmatrix(1,nbasis,1,nbasis);
   e0 = dvector(1,nbasis);
   em = dvector(1,nbasis);
   er = dvector(1,nbasis);
   ae = dvector(1,nbasis);
   be = dvector(1,nbasis);
   ab1 = dvector(1,nbasis);

   clight = 137.036;
   clight2 = clight*clight;
   lang2 = lang*(lang+1.0);
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
      {
      ib = bss[i].ib;                                   //associated bsplines index for i
      jb = bss[j].ib;                                   //associated bsplines index for j
      tmin = MAX(bss[i].tmin,bss[j].tmin);
      tmax = MIN(bss[i].tmax,bss[j].tmax);
      if(tmin >= tmax) pp0[i][j] = ss[i][j] = vv[i][j] = 0.0;
      else
         {
         dt = tmax - tmin;
         sum = sums = sumv = 0.0;
	 if(lang==0) for(q=0;q<QMAX_GAUSS;q++)			//case S-wave: lang = 0.0
            {
            t = (tmin + tmax + dt*egl[q])/2.0;
            dbidbj = dbspline(tn,ib,kb,t)*dbspline(tn,jb,kb,t);
            bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
            sum += wgl[q]*(dbidbj);
            sums += wgl[q]*bibj;
            sumv -= wgl[q]*bibj*zat/t;
            }
         else for(q=0;q<QMAX_GAUSS;q++)
            {
            t = (tmin + tmax + dt*egl[q])/2.0;
            dbidbj = dbspline(tn,ib,kb,t)*dbspline(tn,jb,kb,t);
            bibj = bspline(tn,ib,kb,t)*bspline(tn,jb,kb,t);
            sum += wgl[q]*(dbidbj + bibj*lang2/(t*t));
            sums += wgl[q]*bibj;
            sumv -= wgl[q]*bibj*zat/t;
            }
         pp0[i][j] = sum*dt/2.0;
	 ss[i][j] = sums*dt/2.0;
	 vv[i][j] = sumv*dt/2.0;
         }
      pp0[j][i] = pp0[i][j];					//classical momentum square in Bspline basis: P^2 (T0 = P^2/2
      ss[j][i] = ss[i][j];					//overlap matrix in Bspline basis
      vv[j][i] = vv[i][j];					//potential energy matrix V in Bspline basis
      }


//-------------------------
   diagoz(nbasis,ss,pp0,e0,cc0);				//get P^2 basis: cc0 --- e0[n] = <n|P^2|n>
printf("lang = %d \t zat = %+1.3f \t kb=%d \n",lang,zat,kb);
   for(n=1;n<=nbasis;n++)
      {
      em[n] = sqrt(e0[n]*clight2 + clight2*clight2);		//relativistic free energy
      er[n] = em[n] - clight2;					//relativistic kinetic energy
      ae[n] = sqrt(0.5*(em[n]+clight2)/em[n]);			//A factor (close to 1.0)
      be[n] = clight/(em[n]+clight2);				//B factor (close to 1/2c)
if(n%10 == 0) printf("n=%d \t t0=%+le \t tr=%+le \t ae=%+le \t be=%+le \n",n,e0[n]/2.0,er[n],ae[n],be[n]);
      }

    /*for(n=1;n<=nbasis;n++)  for(m=1;m<=nbasis;m++)
    {
        vvd[n][m] = 0.0;
        for(i=1;i<=nbasis;i++)
        {
            for(j=1;j<=nbasis;j++)
            {
               vvd[n][m] +=  cc0[i][n]*cc0[j][n]*vv[i][j];
            }

        }
        if(n==m) printf("%+le\n",vvd[n][m]);
    }
    exit(1);*/



//-------------------------
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)		//set relativistic kinetic energy in P^2 basis
      {
      kk[i][j] = aa[i][j] = 0.0;
//      for(n=1;n<=nbasis;n++)  kk[i][j] += cc0[i][n]*cc0[j][n]*e0[n]/2.0;	//classical energy
      for(n=1;n<=nbasis;n++)
         {
         kk[i][j] += cc0[i][n]*cc0[j][n]*er[n];			//relativistic energy
         aa[i][j] += cc0[i][n]*cc0[j][n]*ae[n];			//ae in Bspline basis
	 }
      kk[j][i] = kk[i][j];
      aa[j][i] = aa[i][j];
//if(i<=20 && j<=(i+12)) printf("%d  %d\t %+le\n",i,j,ss[i][j]);
      }

    double **ava1, **ava2;
    double **temp,**temp2,**temp3;

    ava1 = dmatrix(1,nbasis,1,nbasis);
    ava2 = dmatrix(1,nbasis,1,nbasis);
    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);
    temp3 = dmatrix(1,nbasis,1,nbasis);



printf("end K\n");
   sstransform(nbasis,kb,ss,kk,ttr);				//back in bsplines basis

//-------------------------
   aatransform(nbasis,aa,vv,avva);				    //A.V.A in Bspline basis
   sstransform(nbasis,kb,ss,avva,vvr);				//back in Bsplines basis
    
    /*for(n=1;n<=nbasis;n++)  for(m=1;m<=nbasis;m++)
    {
        if(n==m) printf("%+le \t %+le \t %+le\n",vvr[n][m],ava1[n][m],ava2[n][m]);
    }

    exit(1);*/

//-------------------------
   if(lang==0)
      {
      darwin = 0.5*zat*dbspline(tn,1,kb,0.0)*dbspline(tn,1,kb,0.0);
      for(i=1;i<=nbasis;i++)
         {
         ab1[i] = 0.0;;
         for(m=1;m<=nbasis;m++) ab1[i] += cc0[i][m]*ae[m]*be[m]*cc0[1][m];
	 }
      for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)		//set Darwin term
         {
         pvvp[i][j] = darwin*ab1[i]*ab1[j];			//Darwin term in P^2 bassis
	 pvvp[j][i] = pvvp[i][j];
	 }
      sstransform(nbasis,kb,ss,pvvp,vvd);			//back in Bsplines basis
      }
   else for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++) vvd[i][j] = 0.0;
     /* {
      for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)		//set Darwin term
         {
         zz[i][j] = yy[i][j] = 0.0;
         for(m=1;m<=nbasis;m++)
            {
            zz[i][j] += cc0[i][m]*ae[m]*be[m]*cc0[j][m];
            yy[i][j] += cc0[i][m]*ae[m]*be[m]*e0[m]*cc0[j][m];
            }
         zz[j][i] = zz[i][j];
         yy[j][i] = yy[i][j];
         }
for(i=1;i<=100;i++) printf("i=%d \t vii=%+le \t zii=%+le \t yii=%+le \t aii=%+le \n",i,vv[i][i],zz[i][i],yy[i][i],aa[i][i]);
      zytransform(nbasis,zz,yy,vv,pvvp);			//(Z.V.Y + Y.V.Z)/2 in Bspline basis
      sstransform(nbasis,kb,ss,pvvp,vvd);			//back in Bsplines basis
      }
double **vvd1;
vvd1 = dmatrix(1,nbasis,1,nbasis);
vvtest(nbasis,kb,ss,cc0,ae,be,e0,vv,vvd1);
for(i=1;i<=150;i++) printf("i=%d \t vdii=%+le \t vd1ii=%+le \t vii=%+le\n",i,vvd[i][i],vvd1[i][i],vv[i][i]);
free_dmatrix(vvd1,1,nbasis,1,nbasis);
getchar();
exit(1);*/

//-------------------------
   for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)
      {
      hh[i][j] = ttr[i][j] + vvr[i][j] + 0.0*vvd[i][j] - 0.5*pp0[i][j] - vv[i][j];	//difference Hamiltonain
//      hh[i][j] = ttr[i][j] - 0.5*pp0[i][j];	//difference Hamiltonian
      hh[j][i] = hh[i][j];
      }

//-------------------------
   if(iprint)
      {
      e1 = dvector(1,nbasis);
      diagoz(nbasis,ss,ttr,e1,cc0);
      for(n=1;n<=MIN(nbasis,100);n++)
         {
//         printf("n=%d \t \t ae=%+le \t t0=%+le \t tr=%+le \n",n,ae[n],e0[n]/2.0,er[n]);
         printf("n=%d \t t0=%+le \t tr=%+le \t vv0=%+le \t vvr=%+le \t vvd=%+le \t hh=%+le\n",n,0.5*pp0[n][n],ttr[n][n],vv[n][n],vvr[n][n],vvd[n][n],hh[n][n]);
//         printf("n=%d\thh0=%+le \t hhr=%+le hh=%+le\n",n,0.5*pp0[n][n]+vv[n][n],ttr[n][n]+vvr[n][n],hh[n][n]);
         }
      free_dvector(e1,1,nbasis);
      }

//---
   free_dmatrix(ss,1,nbasis,1,nbasis);
   free_dmatrix(pp0,1,nbasis,1,nbasis);
   free_dmatrix(ttr,1,nbasis,1,nbasis);
   free_dmatrix(vv,1,nbasis,1,nbasis);
   free_dmatrix(vvr,1,nbasis,1,nbasis);
   free_dmatrix(vvd,1,nbasis,1,nbasis);
   free_dmatrix(avva,1,nbasis,1,nbasis);
   free_dmatrix(pvvp,1,nbasis,1,nbasis);
   free_dmatrix(cc0,1,nbasis,1,nbasis);
   free_dmatrix(aa,1,nbasis,1,nbasis);
   free_dmatrix(kk,1,nbasis,1,nbasis);
   free_dmatrix(ks,1,nbasis,1,nbasis);
   free_dmatrix(yy,1,nbasis,1,nbasis);
   free_dmatrix(zz,1,nbasis,1,nbasis);
   free_dvector(e0,1,nbasis);
   free_dvector(em,1,nbasis);
   free_dvector(er,1,nbasis);
   free_dvector(ae,1,nbasis);
   free_dvector(be,1,nbasis);
   free_dvector(ab1,1,nbasis);
//---
//exit(1);
}



/*-----------------------------------------------------------------
   Transforms a matrix M in P^2 basis to basis spline basis by
   application of L.M.R.
   Matrix ss should be symmetric and band diagonal in [i-kb:i+kb].
-----------------------------------------------------------------*/
void sstransform(int nbasis, int kb, double **ss, double **mm, double **lmmr)
{
    int n,nmin,nmax;;
    int k,kmin,kmax;
    int i,j;
    double **ks;

    ks = dmatrix(1,nbasis,1,nbasis);

    for(n=1;n<=nbasis;n++) for(j=1;j<=nbasis;j++)	//right transform
        {
        ks[n][j] = 0.0;
        kmin = MAX(1,j-kb);
        kmax = MIN(j+kb,nbasis);
        for(k=kmin;k<=kmax;k++) ks[n][j] += mm[n][k]*ss[j][k];
        }
    //printf("end KR\n");
    for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)	//left transform
        {
        lmmr[i][j] = 0.0;
        nmin = MAX(1,i-kb);
        nmax = MIN(i+kb,nbasis);
        for(n=nmin;n<=nmax;n++) lmmr[i][j] += ss[i][n]*ks[n][j];
        lmmr[j][i] = lmmr[i][j];				//tranformed matrix
    //if(j%100==0) printf("i j = %d  %d \n",i,j);
        }
    //printf("end LKR\n");

    free_dmatrix(ks,1,nbasis,1,nbasis);
}



/*-----------------------------------------------------------------
   Transforms a matrix M in P^2 basis to basis spline basis by
   application of A.M.A.
   A should be symmetric.
-----------------------------------------------------------------*/
void aatransform(int nbasis, double **aa, double **mm, double **lmmr)
{
    int n,k;
    int i,j;
    double **ks;

    ks = dmatrix(1,nbasis,1,nbasis);

    for(n=1;n<=nbasis;n++) for(j=1;j<=nbasis;j++)	//right transform
    {
        ks[n][j] = 0.0;
        for(k=1;k<=nbasis;k++) ks[n][j] += mm[n][k]*aa[j][k];
    }
    //printf("end MA\n");
    for(i=1;i<=nbasis;i++) for(j=i;j<=nbasis;j++)	//left transform
    {
        lmmr[i][j] = 0.0;
        for(n=1;n<=nbasis;n++) lmmr[i][j] += aa[i][n]*ks[n][j];
        lmmr[j][i] = lmmr[i][j];				//tranformed matrix
    //if(j%100==0) printf("i j = %d  %d \n",i,j);
    }
    //printf("end AMA\n");

    free_dmatrix(ks,1,nbasis,1,nbasis);
}
