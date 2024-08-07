double c_ua = 137.0;

/* ---------------------
    Base transform tools
------------------------ */
void DKH_transform(int nbasis, double **omega, double **ss, double **mat_to_transf, double **DKH_mat)
{
    int i,j;
    double **omegaT;
    double **temp;

    omegaT = dmatrix(1,nbasis,1,nbasis);
    temp = dmatrix(1,nbasis,1,nbasis);
    
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) omegaT[i][j] = temp[i][j] = 0.0;
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) omegaT[i][j] = omega[j][i];

    mult(mat_to_transf,omega,nbasis,temp);
    mult(omegaT,temp,nbasis,DKH_mat);


    free_dmatrix(omegaT,1,nbasis,1,nbasis);
    free_dmatrix(temp,1,nbasis,1,nbasis);
}

void Standard_transform(int nbasis, double **omega, double **ss, double **mat_to_transf, double **Standard_mat)
{
    int i,j,n;
    double **omegaT;
    double **temp, **temp2, **temp3;

    omegaT = dmatrix(1,nbasis,1,nbasis);
    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);
    temp3 = dmatrix(1,nbasis,1,nbasis);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) omegaT[i][j] = temp[i][j] = temp2[i][j] = temp3[i][j] = 0.0;
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) omegaT[i][j] = omega[j][i];

    mult(omegaT,ss,nbasis,temp);
    mult(mat_to_transf,temp,nbasis,temp2);
    mult(omega,temp2,nbasis,temp3);
    mult(ss,temp3,nbasis,Standard_mat);

    free_dmatrix(omegaT,1,nbasis,1,nbasis);
    free_dmatrix(temp,1,nbasis,1,nbasis);
    free_dmatrix(temp2,1,nbasis,1,nbasis);
    free_dmatrix(temp3,1,nbasis,1,nbasis);
}


/* -----------------------------------------
    Coefficient A/R for DKH transformation
-------------------------------------------- */
void DKH_coeff(int nbasis, double **omega, double **ss, double *eorb, double **Ap, double **Rp)
{
    int i,j;
    double Ep;

    double **temp;

    temp = dmatrix(1,nbasis,1,nbasis);
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) Ap[i][j] = Rp[i][j] = 0.0;
    for(i=1;i<=nbasis;i++)
    {
        Ep = sqrt(2.0*eorb[i]*c_ua*c_ua + c_ua*c_ua*c_ua*c_ua);
        Rp[i][i] = c_ua/(Ep+c_ua*c_ua);
        Ap[i][i] = sqrt((Ep+c_ua*c_ua)/(2.0*Ep));
    }

}


/* ------------------------------------------------
    Calculus of the Kinetic Matrix in the p²-basis 
--------------------------------------------------- */
void DKH_Kinetic_construct(int nbasis, double *eorb, double **dkh_kinetic)
{
    int i,j;

    for(i=1;i<=nbasis;i++)
    {   
        for(j=1;j<=nbasis;j++) dkh_kinetic[i][j] = 0.0;
        dkh_kinetic[i][i] = sqrt(2.0*eorb[i]*c_ua*c_ua + c_ua*c_ua*c_ua*c_ua);
    }
}

/*  --------------------------------------------------------------------------------------------------------------------------------------------
    Calculus of the Ionique Matrix in the p²-basis. 
    Ionique_calculus and Ionique_potential work together. Ionique_potential calls Ionique_calculus in order to treat with the C2 coefficient.
    Ionique_construct construct the Ionique Matrix in the p²-basis, in agreements with the quasi-relativistic DKH method.  
    Ionique'{p²} = A * (Ionique + R * Ionique * R) * A                     
-----------------------------------------------------------------------------------------------------------------------------------------------  */

void DKH_Ionique_construct(int nbasis, double **Ap, double **Rp, double **ionique, double **coulomb, double **pVp, double **dkh_ionique, double **omega, double **ss)
{
    int i,j;
    double **temp, **temp2, **temp3, **ctemp, **ctemp2;

    ctemp = dmatrix(1,nbasis,1,nbasis);
    ctemp2 = dmatrix(1,nbasis,1,nbasis);

    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);
    temp3 = dmatrix(1,nbasis,1,nbasis);
    
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp[i][j] = temp2[i][j] = temp3[i][j] = ctemp[i][j] = ctemp2[i][j] = 0.0;
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp2[i][j] = ionique[i][j] + coulomb[i][j];
    DKH_transform(nbasis,omega,ss,temp2,temp); 
    DKH_transform(nbasis,omega,ss,pVp,ctemp2);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++)  temp2[i][j] = 0.0;
    mult(ctemp2,Rp,nbasis,temp2);
    mult(Rp,temp2,nbasis,temp3);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp2[i][j] = 0.0;
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp2[i][j] = temp[i][j] + temp3[i][j];  //The minus come from the double operator p = - i d/dr

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp[i][j] =  ctemp[i][j] = dkh_ionique[i][j] = 0.0;
    mult(temp2,Ap,nbasis,temp);
    mult(Ap,temp,nbasis,dkh_ionique);

    free_dmatrix(ctemp,1,nbasis,1,nbasis);
    free_dmatrix(ctemp2,1,nbasis,1,nbasis);
    free_dmatrix(temp,1,nbasis,1,nbasis);
    free_dmatrix(temp2,1,nbasis,1,nbasis);
    free_dmatrix(temp3,1,nbasis,1,nbasis);
}

/*  --------------------------------------------------------------------------------------------------------------------------------------------
    Calculus of the Odd-Matrix (off-diagonal matrix) in the p²-basis. 
    O1_calculus and O1_potential work together. O1_potential calls O1_calculus in order to treat with the C2 coefficient.
    O1_construct construct the off-diagonal matrix in the p²-basis, in agreements with the quasi-relativistic DKH method.
    O1{p²} = A * (R * Ionique - Ionique * R) * A                      
-----------------------------------------------------------------------------------------------------------------------------------------------  */
void DKH_O1_construct(int nbasis, double **Ap, double **Rp, double **pV, double **Vp, double **dkh_odd, double **omega, double **ss)
{
    int i,j,k;
    double **temp, **temp2, **mtemp, **mtemp2;
    double **ttemp, **ttemp2;
    
    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);
    ttemp = dmatrix(1,nbasis,1,nbasis);
    ttemp2 = dmatrix(1,nbasis,1,nbasis);

    mtemp = dmatrix(1,nbasis,1,nbasis);
    mtemp2 = dmatrix(1,nbasis,1,nbasis);
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp[i][j] = temp2[i][j] = mtemp[i][j] = mtemp2[i][j] = 0.0;

    DKH_transform(nbasis,omega,ss,pV,temp);
    DKH_transform(nbasis,omega,ss,Vp,temp2);

    mult(Rp,pV,nbasis,mtemp);
    mult(Vp,Rp,nbasis,mtemp2);
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) ttemp[i][j] = mtemp[i][j] - mtemp2[i][j];

    mult(Ap,ttemp,nbasis,ttemp2);
    mult(ttemp2,Ap,nbasis,dkh_odd);

    free_dmatrix(temp,1,nbasis,1,nbasis);
    free_dmatrix(ttemp,1,nbasis,1,nbasis);
    free_dmatrix(temp2,1,nbasis,1,nbasis);
    free_dmatrix(ttemp2,1,nbasis,1,nbasis);
    free_dmatrix(mtemp,1,nbasis,1,nbasis);
    free_dmatrix(mtemp2,1,nbasis,1,nbasis);
}


/*  --------------------------------------------------------------------------------------------------------------------------------------------
    Calculus of the W-Matrix (off-diagonal matrix) in the p²-basis, in agreements with the quasi-relativistic DKH method. 
    W1[i][j] = O1[i][j] / ( E[i] + E[j] )                      
-----------------------------------------------------------------------------------------------------------------------------------------------  */
void DKH_Kernel_W(int nbasis, double *eorb, double **omega, double **ss, double **dkh_odd, double **dkh_kernel)
{
    int i,j;
    double Ei, Ej;
    double **temp, **temp2;

    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp[i][j] = temp2[i][j] = 0.0;

    DKH_transform(nbasis,omega,ss,dkh_odd,temp);

    for(i=1;i<=nbasis;i++)
    {
        Ei = sqrt(2.0*eorb[i]*c_ua*c_ua + c_ua*c_ua*c_ua*c_ua);
        for(j=1;j<=nbasis;j++)
        {
            Ej = sqrt(2.0*eorb[j]*c_ua*c_ua + c_ua*c_ua*c_ua*c_ua);
            //temp2[i][j] = temp[i][j]/(Ei+Ej);
            dkh_kernel[i][j] = temp[i][j]/(Ei+Ej);
        }
    }

    //Standard_transform(nbasis,omega,ss,temp2,dkh_kernel);

    free_dmatrix(temp,1,nbasis,1,nbasis);
    free_dmatrix(temp2,1,nbasis,1,nbasis);
}



/* ------------------------------------------------
    Construction of the 2nd order of DKH
--------------------------------------------------- */
void DKH_HDK2(int nbasis, double **dkh_eshift, double **dkh_kinetic, double **dkh_ionique, double **dkh_odd, double **dkh_kernel, double **dkh_hdk2, double **omega, double **ss)
{
    int i,j,k;
    double **temp, **temp2, **temp3;
    double **mtemp, **mtemp2;
    double **HDK1;


    temp = dmatrix(1,nbasis,1,nbasis);
    temp2 = dmatrix(1,nbasis,1,nbasis);
    temp3 = dmatrix(1,nbasis,1,nbasis);
    mtemp = dmatrix(1,nbasis,1,nbasis);
    mtemp2 = dmatrix(1,nbasis,1,nbasis);
    HDK1 = dmatrix(1,nbasis,1,nbasis);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp[i][j] = temp2[i][j] = temp3[i][j] = mtemp[i][j] = mtemp2[i][j] = HDK1[i][j] = 0.0;
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) HDK1[i][j] = dkh_kinetic[i][j] + dkh_ionique[i][j] + dkh_eshift[i][j];
    //Standard_transform(nbasis,omega,ss,HDK1,temp2);

    mult(dkh_kernel,dkh_odd,nbasis,mtemp);
    mult(dkh_odd,dkh_kernel,nbasis,mtemp2);

    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) temp3[i][j] = -(mtemp[i][j] + mtemp2[i][j]);  //minus come from the i in p = - i d/dr when we compute O-matrix
    for(i=1;i<=nbasis;i++) for(j=1;j<=nbasis;j++) dkh_hdk2[i][j] = HDK1[i][j] - 0.5*temp3[i][j];

    free_dmatrix(temp,1,nbasis,1,nbasis);
    free_dmatrix(temp2,1,nbasis,1,nbasis);
    free_dmatrix(temp3,1,nbasis,1,nbasis);

    free_dmatrix(mtemp,1,nbasis,1,nbasis);
    free_dmatrix(mtemp2,1,nbasis,1,nbasis);
    free_dmatrix(HDK1,1,nbasis,1,nbasis);
}


void DKH_construct(int nbasis, int kb, double *tn, int lang, double **ss, double **kinetic, double **ionique, double **coulomb, double **dkh_hdk2, struct bsp_basis *bss, struct vxaimp *pspa)
{
    int i,j,n;

    double *eorb, **corb;
    double **omega;

    double **Ap, **Rp;
    double **dkh_kinetic, **dkh_ionique, **dkh_odd, **dkh_kernel;
    double **dkh_eshift;
    double **dk2temp, **dk3temp;

    double **pv, **vp, **pvp;
    double **kinetic_test;

    eorb = dvector(1,nbasis);							
    corb = dmatrix(1,nbasis,1,nbasis);
    omega = dmatrix(1,nbasis,1,nbasis);

    Ap = dmatrix(1,nbasis,1,nbasis);
    Rp = dmatrix(1,nbasis,1,nbasis);

    dkh_kinetic = dmatrix(1,nbasis,1,nbasis);
    dkh_ionique = dmatrix(1,nbasis,1,nbasis);
    dkh_odd = dmatrix(1,nbasis,1,nbasis);
    dkh_kernel = dmatrix(1,nbasis,1,nbasis);

    dkh_eshift = dmatrix(1,nbasis,1,nbasis);
    kinetic_test = dmatrix(1,nbasis,1,nbasis);

    dk2temp = dmatrix(1,nbasis,1,nbasis);
    dk3temp = dmatrix(1,nbasis,1,nbasis);

    pv = dmatrix(1,nbasis,1,nbasis);
    vp = dmatrix(1,nbasis,1,nbasis);
    pvp = dmatrix(1,nbasis,1,nbasis);


    for(i=1;i<=nbasis;i++)
    {   
        eorb[i] = 0.0;
        for(n=1;n<=nbasis;n++) corb[i][n] = omega[i][n] = dkh_eshift[i][n] = 0.0;
        dkh_eshift[i][i] -= c_ua*c_ua;
    }

    diagoz(nbasis,ss,kinetic,eorb,corb);
    for(i=1;i<=nbasis;i++) for(n=1;n<=nbasis;n++) omega[i][n] = corb[i][n];


    DKH_coeff(nbasis,omega,ss,eorb,Ap,Rp);
    DKH_Kinetic_construct(nbasis,eorb,dkh_kinetic);


    dmkvvl(nbasis,kb,tn,lang,bss,pspa,pv,vp,pvp);
    DKH_Ionique_construct(nbasis,Ap,Rp,ionique,coulomb,pvp,dkh_ionique,omega,ss);

    DKH_O1_construct(nbasis,Ap,Rp,pv,vp,dkh_odd,omega,ss);
    DKH_Kernel_W(nbasis,eorb,omega,ss,dkh_odd,dkh_kernel);
    DKH_HDK2(nbasis,dkh_eshift,dkh_kinetic,dkh_ionique,dkh_odd,dkh_kernel,dk2temp,omega,ss);
    Standard_transform(nbasis,omega,ss,dk2temp,dkh_hdk2);
    
    /*for(i=1;i<=nbasis;i++) printf("%+le\n",dkh_hdk2[i][i]);
    printf("\n\n\n");*/


}
