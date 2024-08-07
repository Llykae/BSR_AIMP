int found_occu(char *occu, FILE *fin); 
void init(struct gauss_basis_set *gbs,  char *data, char *atom);
void readgbs(struct gauss_basis_set *gbs, struct orbital *orb, char *data, char *atom, int swc);

double gtonorm(double a, int px, int py, int pz);
double gtoovrlp(struct orbital *bss, int i, int j, int k, int h);


/* Return the line of the first occurence of a word  */ 

int found_occu(char *occu, FILE *fin)     
{
  char chaine[200] = "";
  int nLines = 1;

  if (fin != NULL){
    while (fgets(chaine, 200, fin) != NULL)
      {
	if (strstr(chaine, occu) != NULL){return nLines;}
	else{nLines += 1;}}}
  else{printf("ERROR FILE");}

  return nLines;
}


/* Read the size of the basis set :
   ---> ne, number of electrons
   ---> ngs, number of s-type gaussian
   ---> nfs, number of s-type orbital   */

void init(struct gauss_basis_set *gbs,  char *data, char *atom)
{
  int i,swc;
  int nLine = 0;
  char trash[255];
  FILE *fin;

  gbs[0].nfs = 0;
  gbs[0].nfp = 0;
  gbs[0].nfd = 0;
  
  fin = fopen(data, "r");
  nLine = found_occu(atom, fin);
  rewind(fin);
  
  for(i=0;i<nLine;i++){fgets(trash, 255, fin);}
  fscanf(fin,"%lf",&gbs[0].ne);                  
  rewind(fin);

  if((gbs[0].ne > 4.0) && (gbs[0].ne <= 18.0))      {swc = 1;}   //p -> higher core orbital
  else if(gbs[0].ne > 18.0)                         {swc = 2;}   //d -> higher core orbital
  
  nLine = found_occu("S-type",fin);
  rewind(fin);

  for(i=0;i<nLine;i++){fgets(trash,255,fin);}
  fscanf(fin,"%d %d",&gbs[0].ngs,&gbs[0].nfs);
  rewind(fin);
  switch(swc){
  case 1 :

    nLine = found_occu("P-type",fin);
    rewind(fin);

    for(i=0;i<nLine;i++){fgets(trash,255,fin);}
    fscanf(fin,"%d %d",&gbs[0].ngp,&gbs[0].nfp);
    
    rewind(fin);
    break;

  case 2 :

    nLine = found_occu("P-type",fin);
    rewind(fin);

    for(i=0;i<nLine;i++){fgets(trash,255,fin);}
    fscanf(fin,"%d %d",&gbs[0].ngp,&gbs[0].nfp);
    rewind(fin);

    nLine = found_occu("D-type",fin);
    rewind(fin);

    for(i=0;i<nLine;i++){fgets(trash,255,fin);}
    fscanf(fin,"%d %d",&gbs[0].ngd,&gbs[0].nfd);
    rewind(fin);
    break;}

  fclose(fin);
}

/* Read the Gaussian basis set of the core-orbital */

void readgbs(struct gauss_basis_set *gbs, struct orbital *orb, char *data, char *atom, int swc)
{
  int i,j,k,count;
  int nLine = 0;
  double sums[gbs[0].nfs],sump[gbs[0].nfp],sumd[gbs[0].nfd];
  
  char trash[255], t[255];
  char *t2 = NULL;
  FILE *fin;

  fin = fopen(data, "r");
  nLine = found_occu("S-type",fin);
  rewind(fin);

  /*-------- nS ORBITAL BASIS SET --------*/
  orb[1].l = 0;
  for(i=0;i<nLine+1;i++){fgets(trash,255,fin);}
  for(i=1;i<=gbs[0].ngs;i++){fscanf(fin,"%lf", &(orb[1].arg[i]));} 		
  for(i=1;i<=gbs[0].ngs;i++){orb[1].nrm[i] = gtonorm(orb[1].arg[i],0,0,0);} 		
  for(j=2;j<=gbs[0].nfs;j++) for(i=1;i<=gbs[0].ngs;i++){
	orb[j].arg[i] = orb[1].arg[i];
	orb[j].nrm[i] = orb[1].nrm[i];
	orb[j].l = orb[1].l;}

  for(i=0;i<=gbs[0].ngs;i++){
    fgets(t,255,fin);
    t2 = strtok(t," ");
    j = 1;
    while(t2 != NULL){
      sscanf(t2,"%lf",&(orb[j].cgauss[i]));
      j += 1;
      t2 = strtok(NULL," "); 
    }     
  }

  for(k=1;k<=gbs[0].nfs;k++){
    sums[k]=0.0;
    for(i=1;i<=gbs[0].ngs;i++){
      for(j=1;j<=gbs[0].ngs;j++){sums[k]+=orb[k].cgauss[i]*orb[k].cgauss[j]*gtoovrlp(orb,i,j,k,k);}}}

  for(k=1;k<=gbs[0].nfs;k++){
    for(i=1;i<=gbs[0].ngs;i++){
      orb[k].cgauss[i] /= sqrt(sums[k]);}}
  
  rewind(fin);
  
  switch(swc){
  case 1 :

    /*------- nP ORBITAL BS --------*/
    nLine = found_occu("P-type",fin);
    rewind(fin);

    count = gbs[0].nfs;
    orb[1+count].l = 1;
    for(i=0;i<nLine+1;i++){fgets(trash,255,fin);}
    for(i=1;i<=gbs[0].ngp;i++){fscanf(fin,"%lf", &(orb[1+count].arg[i]));} 
    for(i=1;i<=gbs[0].ngp;i++){orb[1+count].nrm[i] = gtonorm(orb[1+count].arg[i],1,0,0);}
    for(j=2;j<=gbs[0].nfp;j++) for(i=1;i<=gbs[0].ngp;i++){
        orb[j+count].arg[i] = orb[1+count].arg[i];
        orb[j+count].nrm[i] = orb[1+count].nrm[i];
	orb[j+count].l = 1;}


    for(i=0;i<=gbs[0].ngp;i++){
      fgets(t,255,fin);
      t2 = strtok(t," ");
      j = 1;
      while(t2 != NULL){
	sscanf(t2,"%lf",&(orb[j+count].cgauss[i]));
	j += 1;
	t2 = strtok(NULL," "); 
      }
    }
    rewind(fin);

    for(k=1;k<=gbs[0].nfp;k++){
      sump[k]=0.0;
      for(i=1;i<=gbs[0].ngp;i++){
	for(j=1;j<=gbs[0].ngp;j++){
	  sump[k]+=orb[k+count].cgauss[i]*orb[k+count].cgauss[j]*gtoovrlp(orb,i,j,k+count,k+count);}}}

    for(k=1;k<=gbs[0].nfp;k++){
      for(i=1;i<=gbs[0].ngp;i++){
	orb[k+count].cgauss[i] /= sqrt(sump[k]);}}
    break;

  case 2 :

    /*------- nP ORBITAL BS --------*/
    nLine = found_occu("P-type",fin);
    rewind(fin);

    count = gbs[0].nfs;
    orb[1+count].l = 1;
    for(i=0;i<nLine+1;i++){fgets(trash,255,fin);}
    for(i=1;i<=gbs[0].ngp;i++){fscanf(fin,"%lf", &(orb[1+count].arg[i]));} 
      for(i=1;i<=gbs[0].ngp;i++){orb[1+count].nrm[i] = gtonorm(orb[1+count].arg[i],1,0,0);}
      for(j=2;j<=gbs[0].nfp;j++) for(i=1;i<=gbs[0].ngp;i++){
	  orb[j+count].arg[i] = orb[1+count].arg[i];
	  orb[j+count].nrm[i] = orb[1+count].nrm[i];
	  orb[j+count].l = orb[1+count].l;}


      for(i=0;i<=gbs[0].ngp;i++){
	fgets(t,255,fin);
	t2 = strtok(t," ");
	j = 1;
	while(t2 != NULL){
	  sscanf(t2,"%lf",&(orb[j+count].cgauss[i]));
	  j += 1;
	  t2 = strtok(NULL," "); 
	}
      }

      for(k=1;k<=gbs[0].nfp;k++){
	sump[k]=0.0;
	for(i=1;i<=gbs[0].ngp;i++){
	  for(j=1;j<=gbs[0].ngp;j++){
	    sump[k]+=orb[k+count].cgauss[i]*orb[k+count].cgauss[j]*gtoovrlp(orb,i,j,k+count,k+count);}}}

      for(k=1;k<=gbs[0].nfp;k++){
	for(i=1;i<=gbs[0].ngp;i++){
	  orb[k+count].cgauss[i] /= sqrt(sump[k]);}}
      rewind(fin);
    
      /*------- nD ORBITAL BS --------*/
      nLine = found_occu("D-type",fin);
      rewind(fin);

      count += gbs[0].nfp;
      orb[1+count].l = 2;
      for(i=0;i<nLine+1;i++){fgets(trash,255,fin);}
      for(i=1;i<=gbs[0].ngd;i++){fscanf(fin,"%lf", &(orb[1+count].arg[i]));}
      for(i=1;i<=gbs[0].ngd;i++){orb[1+count].nrm[i] = gtonorm(orb[1+count].arg[i],2,0,0);}
      for(j=2;j<=gbs[0].nfd;j++) for(i=1;i<=gbs[0].ngd;i++){
	  orb[j+count].arg[i] = orb[1+count].arg[i];
	  orb[j+count].nrm[i] = orb[1+count].nrm[i];
	  orb[j+count].l = orb[1+count].l;}

      for(i=0;i<=gbs[0].ngd;i++){
	fgets(t,255,fin);
	t2 = strtok(t," ");
	j = 1;
	while(t2 != NULL){
	  sscanf(t2,"%lf",&(orb[j+count].cgauss[i]));
	  j += 1;
	  t2 = strtok(NULL," "); 
	}
      }


      for(k=1;k<=gbs[0].nfd;k++){
	sumd[k]=0.0;
	for(i=1;i<=gbs[0].ngd;i++){
	  for(j=1;j<=gbs[0].ngd;j++){
	    sumd[k]+=orb[k+count].cgauss[i]*orb[k+count].cgauss[j]*gtoovrlp(orb,i,j,k+count,k+count);}}}

      for(k=1;k<=gbs[0].nfd;k++){
	for(i=1;i<=gbs[0].ngd;i++){
	  orb[k+count].cgauss[i] /= sqrt(sumd[k]);}}
      rewind(fin);
      break;
  }

  fclose(fin);
}

