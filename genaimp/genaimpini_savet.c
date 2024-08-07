/*-------------------------------------------------------------------------------
------------------------------------------------------------------------------*/


double gtonorm(double a, int p1, int p2, int p3);

void initdata(int *natm, int *nelc, int *nsym);
void initatm(int ntot, struct atom_dat **at_dt);
void initmath();



/*-----------------------------------------------------------------------------
   initdata() reads:
        -the number of atoms in the cluster: natm
        -the number of electrons in the cluster: nelc
        -the number of symmetries: nsym
------------------------------------------------------------------------------*/
void initdata(int *natm, int *nelc, int *nsym)
{
   FILE *fdin;
   char strtmp[250];

   fdin=fopen("input/data.dat","r");
   fgets(strtmp,200,fdin);
   fgets(strtmp,200,fdin);
   fgets(strtmp,200,fdin);
   fscanf(fdin,"%d %s\n",natm,strtmp);
   fscanf(fdin,"%d %s\n",nelc,strtmp);
   fscanf(fdin,"%d %s\n",nsym,strtmp);
   fclose(fdin);
}



/*-----------------------------------------------------------------------------
   initatm() initializes atomic data structure for atom.
-----------------------------------------------------------------------------*/
void initatm(int ntot, struct atom_dat **at_dt)
{
   FILE *fdin;
 char strtmp[256];
 char output[256];
 char strips[16];
 int iat,i;

   (*at_dt) = malloc((ntot+1)*sizeof(struct atom_dat));		//memory for atom data structure at_dt

   fdin=fopen("input/system.dat","r");				//read input file "system.dat"
   for (i=0;i<8;i++) fgets(strtmp,200,fdin);

   for(iat=1;iat<=ntot;iat++)
      {
       fscanf(fdin,"%s %d %d %lf",(*at_dt)[iat].element,&(*at_dt)[iat].imol,&(*at_dt)[iat].zat ,&(*at_dt)[iat].m);				fgets(strtmp,200,fdin);
       fscanf(fdin,"%lf %lf %lf",&(*at_dt)[iat].r[1],&(*at_dt)[iat].r[2],&(*at_dt)[iat].r[3]);							fgets(strtmp,200,fdin);
       fscanf(fdin,"%lf %lf %lf",&(*at_dt)[iat].p[1],&(*at_dt)[iat].p[2],&(*at_dt)[iat].p[3]);							fgets(strtmp,200,fdin);
       fscanf(fdin,"%s",(*at_dt)[iat].bsnm);													fgets(strtmp,200,fdin);
       fscanf(fdin,"%s",(*at_dt)[iat].pspnm);													fgets(strtmp,200,fdin);
// 	 printf("%s %d %d %f\n",	(*at_dt)[iat].element,	(*at_dt)[iat].imol,	(*at_dt)[iat].zat ,	(*at_dt)[iat].m);
// 	 printf("%f\t%f\t%f\n",		(*at_dt)[iat].r[1],	(*at_dt)[iat].r[2],	(*at_dt)[iat].r[3]);
// 	 printf("%f\t%f\t%f\n",		(*at_dt)[iat].p[1],	(*at_dt)[iat].p[2],	(*at_dt)[iat].p[3]);
      printf("%s\n",(*at_dt)[iat].bsnm);
      printf("%s\n",(*at_dt)[iat].pspnm);

      norm((*at_dt)[iat].p,3);
      norm((*at_dt)[iat].r,3);
      }
   fclose(fdin);   
}



/*-----------------------------------------------------------------------------
   initmath() intializes all mathematical constant needed for programm lcao.
-----------------------------------------------------------------------------*/
void initmath()
{
   Pi = pi = acos(-1.0);
   sqpi = sqrt(pi);

   dawson(1.0);					//initialization of static coefficients for dawson()
   intbin();					//intializes the binomial coefficients Cnp
}
