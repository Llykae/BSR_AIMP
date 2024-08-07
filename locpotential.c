/*---------------------------------------------------------------------
        This file contains the derivation of all local potentials
---------------------------------------------------------------------*/

void density(struct gauss_basis_set *gbs, struct orbital *orb, struct atom *at);
double vpotloc(int z, double r, struct gauss_basis_set *gbs, struct orbital *orb);

void density(struct gauss_basis_set *gbs, struct orbital *orb, struct atom *at)
{
  int i,j,k,count;
  double sum=0.0;

  if(gbs[0].nfs>0){
    for(k=1;k<=gbs[0].nfs;k++){
      for(i=1;i<=gbs[0].ngs;i++){
	for(j=i;j<=gbs[0].ngs;j++){
	  orb[k].cdensity[i][j] = orb[k].cgauss[i]*orb[k].cgauss[j];
	  orb[k].cdensity[j][i] = orb[k].cdensity[i][j];}}}

  /*for(i=1;i<=gbs[0].ngs;i++) for(j=1;j<=gbs[0].ngs;j++) ex[i][j]=0.0;
  for(k=1;k<=gbs[0].nfs;k++) for(i=1;i<=gbs[0].ngs;i++) for(j=1;j<=gbs[0].ngs;j++) ex[i][j]+=orb[k].cdensity[i][j];
  for(i=1;i<=gbs[0].ngs;i++){
    for(j=1;j<=gbs[0].ngs;j++){printf("%+le   ", ex[i][j]);}
    printf("\n");}
  printf("\n\n");*/

    for(k=1;k<=gbs[0].nfs;k++){
      for(i=1;i<=gbs[0].ngs;i++){
	for(j=1;j<=gbs[0].ngs;j++){sum+=2.0*orb[k].cdensity[i][j]*gtoovrlp(orb,i,j,k,k);}}}
    at[0].occ[1] = sum;
    printf("Density of S-states : %f\n",sum);}

  if(gbs[0].nfp>0){
    sum = 0.0;
    count = gbs[0].nfs;
    
    for(k=1;k<=gbs[0].nfp;k++){
      for(i=1;i<=gbs[0].ngp;i++){
	for(j=i;j<=gbs[0].ngp;j++){
	  orb[k+count].cdensity[i][j] = orb[k+count].cgauss[i]*orb[k+count].cgauss[j];
	  orb[k+count].cdensity[j][i] = orb[k+count].cdensity[i][j];}}}

  /*ex2=dmatrix(1,gbs[0].ngp,1,gbs[0].ngp);
  for(i=1;i<=gbs[0].ngp;i++) for(j=1;j<=gbs[0].ngp;j++) ex2[i][j]=0.0;
  for(k=1;k<=gbs[0].nfp;k++) for(i=1;i<=gbs[0].ngp;i++) for(j=1;j<=gbs[0].ngp;j++) ex2[i][j]+=orb[k+count].cdensity[i][j];
  for(i=1;i<=gbs[0].ngp;i++){
    for(j=1;j<=gbs[0].ngp;j++){printf("%+le   ", ex2[i][j]);}
    printf("\n");}
  printf("\n\n");

  for(i=1;i<=gbs[0].ngp;i++){
    for(j=1;j<=gbs[0].ngp;j++){printf("%+le   ", 2*ex2[i][j]);}
    printf("\n");}
  exit(1);*/

    for(k=1;k<=gbs[0].nfp;k++){
      for(i=1;i<=gbs[0].ngp;i++){
	for(j=1;j<=gbs[0].ngp;j++){sum+=6.0*orb[k+count].cdensity[i][j]*gtoovrlp(orb,i,j,k+count,k+count);}}}
    at[0].occ[2] = sum;
    printf("Density of P-states : %f\n",sum);}

 
  if(gbs[0].nfd>0){

    sum = 0.0;
    count += gbs[0].nfp;
    
    for(k=1;k<=gbs[0].nfd;k++){
      for(i=1;i<=gbs[0].ngd;i++){
	for(j=i;j<=gbs[0].ngd;j++){
	  orb[k+count].cdensity[i][j] = orb[k+count].cgauss[i]*orb[k+count].cgauss[j];
	  orb[k+count].cdensity[j][i] = orb[k+count].cdensity[i][j];}}}

    for(k=1;k<=gbs[0].nfd;k++){
      for(i=1;i<=gbs[0].ngd;i++){
	for(j=1;j<=gbs[0].ngd;j++){sum+=10.0*orb[k+count].cdensity[i][j]*gtoovrlp(orb,i,j,k+count,k+count);}}}
    at[0].occ[3] = sum;
    printf("Density of D-states : %f\n",sum);}
}

double vpotloc(int z, double r, struct gauss_basis_set *gbs, struct orbital *orb)
{
  int i, j, k, count=0.0;
  double a, mps;
  double vnuc, vs, vp, vd, vloc, vkl=0.0;
  //FILE *fout;

  //fout=fopen("a.out","a");
  
  vnuc = -z/r;
  vs = 0.0;

  for(k=1;k<=gbs[0].nfs;k++){
    for(i=1;i<=gbs[0].ngs;i++){
      for(j=1;j<=gbs[0].ngs;j++){
	a = orb[k].arg[i] + orb[k].arg[j];
	vkl = orb[k].cdensity[i][j]*orb[k].nrm[j]*orb[k].nrm[i]*0.25*sqrt(pi/(a*a*a))/r*erf(sqrt(a)*r);
	vs += vkl;}}}
  vs *= 4*pi*2.0;

  
  count = gbs[0].nfs;
  vp = 0.0;

  for(k=1;k<=gbs[0].nfp;k++){
    for(i=1;i<=gbs[0].ngp;i++){
      for(j=1;j<=gbs[0].ngp;j++){
	a = orb[k+count].arg[i] + orb[k+count].arg[j];
	mps = (3./8)*sqrt(pi/(a*a*a*a*a))*erf(sqrt(a)*r)/r - 1./(4*a*a)*exp(-a*r*r);  //multipolar expansion for l = 0
	vkl = orb[k+count].cdensity[i][j]*orb[k+count].nrm[j]*orb[k+count].nrm[i]*mps;
	vp += vkl;}}}
  vp *= 4*pi*6.0/3.;

  
  count += gbs[0].nfp;
  vd = 0.0;
  mps = 0.0;
  
  for(k=1;k<=gbs[0].nfd;k++){
    for(i=1;i<=gbs[0].ngd;i++){
      for(j=1;j<=gbs[0].ngd;j++){
	a = orb[k+count].arg[i] + orb[k+count].arg[j];
	mps = (15./16)*sqrt(pi/(a*a*a*a*a*a*a))*erf(sqrt(a)*r)/r - 7./(8*a*a*a)*exp(-a*r*r) - 0.25*r*r/(a*a)*exp(-a*r*r);  //multipolar expansion for l = 0
	vkl = orb[k+count].cdensity[i][j]*orb[k+count].nrm[j]*orb[k+count].nrm[i]*mps;
	vd += vkl;}}}
  vd *= 4*pi*10.0/5.;

  vloc =  vnuc + vs + vp + vd;
  return vloc;
}
