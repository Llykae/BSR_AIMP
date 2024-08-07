/*--------------------------------------------------------------------------------

    Construct the matrix G[i][j](L) for the propagation of the wavefunction
    phi_L(t+dt)[i] = exp{-I*H0[i][j](L)*dt} phi_L(t)[j] = G[i][j][L] phi_L(t)[j]
    with H0(L) is the time-independent hamiltonian for the orbiatl L wavefunction 
    
    INPUT:	Lmax = angular momentum  L = 1...Lmax
    		Nmax = dimention or radial grid
    		Vps = time-independent effectif potentiel
    		xt = (N-1) collocation points 
    		
    OUTPUT:     G[i][j][L] propagation matrix for a given L space
    
----------------------------------------------------------------------------------*/
void tred2(double **a, int n, double d[], double e[]) 
{
   int l,k,j,i;
   double scale,hh,h,g,f;
		    	  
   for(i=n;i>=2;i--)
      {
      l=i-1;
      h=scale=0.0;
      if (l>1) 
         {
         for (k=1;k<=l;k++) scale += fabs(a[i][k]);
         if (scale == 0.0) e[i]=a[i][l];
         else 
            {
            for(k=1;k<=l;k++)
               {
               a[i][k] /= scale;
               h += a[i][k]*a[i][k];
               }
            f=a[i][l];
            g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i][l]=f-g;
            f=0.0;
            for(j=1;j<=l;j++)
               {
               a[j][i]=a[i][j]/h;
               g=0.0;
               for(k=1;k<=j;k++)
               g += a[j][k]*a[i][k];
               for(k=j+1;k<=l;k++)
               g += a[k][j]*a[i][k];
               e[j]=g/h;
               f += e[j]*a[i][j];
               }
            hh=f/(h+h);
            for(j=1;j<=l;j++)
               {
               f=a[i][j];
               e[j]=g=e[j]-hh*f;
               for(k=1;k<=j;k++)
               a[j][k] -= (f*e[k]+g*a[i][k]);
               }
            }
         }
      else
      e[i]=a[i][l];
      d[i]=h;
      }
			
   d[1]=0.0;
   e[1]=0.0;
   for(i=1;i<=n;i++)
      {
      l=i-1;
      if(d[i])
         {
         for(j=1;j<=l;j++)
            {
            g=0.0;
            for(k=1;k<=l;k++)
            g += a[i][k]*a[k][j];
            for (k=1;k<=l;k++)
            a[k][j] -= g*a[k][i];
            }
         }
      d[i]=a[i][i];
      a[i][i]=1.0;
      for(j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
      }
}


 
/*-------------------------------------------------------------------------*/ 
void tqli(double d[], double e[], int nn, double **q)
{
   double pythag();
   int i,k,l,m,itr;
   double b,f,g,r,s,c,p,dd;

   for (i=2;i<=nn;i++) e[i-1]=e[i];
   e[nn]=0.0;
   
   for(l=1;l<=nn;l++)
      {
      itr=0;
      do
         {
         for(m=l;m<=nn-1;m++)				/*1-split the matrix for small*/
            {						/*subdiagonal element e.2-tests*/
            dd = fabs(d[m]) + fabs(d[m+1]);		/*on the convergence*/
            if((double)(fabs(e[m]) + dd) == dd) break;
            }
         if(m != l)
            {
            if (itr++ == 30) nerror("too many iterations in tqli");					/*iteration count*/
            g = (d[l+1]-d[l])/(2.0*e[l]);
/*---r = sqrt(1.0 + g*g);---*/
            r = pythag(1.0,g);
            g = d[m] - d[l] + e[l]/(g + SIGN(r,g));   /*g = d[m] - ks*/
            s = c = 1;
            p = 0.0;
            for(i=m-1;i>=l;i--)
               {
               f = s*e[i];
               b = c*e[i];
/*---e[i+1] = (r=sqrt(f*f+g*g));---*/
               e[i+1] = (r=pythag(f,g));
               if(r==0.0)
                  {
                  d[i+1] -= p;
                  e[m] = 0.0;
                  break;
                  }
               s = f/r;
               c = g/r;
               g = d[i+1] - p;
               r = (d[i] - g)*s + 2.0*c*b;
               d[i+1] = g + (p=s*r);
               g = c*r - b;
               for(k=1;k<=nn;k++)                          /*calculate the change*/
                  {                                       /*of frame matrix*/
                  f = q[k][i+1];
                  q[k][i+1] = s*q[k][i] + c*f;
                  q[k][i] = c*q[k][i] - s*f;
                  }
               }
			   
            if(r==0.0 && i>=l) continue;
            d[l]-=p;
            e[l]=g;
            e[m]=0.0;
            }
         } while (m != l);
      }
}



/*-------------------------eigsrt()----------------------------
   eigsrt sorts the eigenvalues into ascending order, and 
   rearranges the columns of v correspondingly.
   See Numerical Recipies in C.
	input :	d[1...n]   the eigenvalues
	        v[1...n][1...n]   the eigenvectors
	
	output : n number of element to sort
-------------------------------------------------------------*/
void eigsrt(double *d, double **v, int n)
{
   int k,j,i;
   double p;
   
   for(i=1;i<n;i++)
      {
      p = d[k=i];
      for(j=i+1;j<=n;j++)
         if(d[j] <= p) p = d[k=j];
         /*if(d[j] >= p) p = d[k=j];*/
      if(k != i)
         {
         d[k] = d[i];
         d[i] = p;
         for(j=1;j<=n;j++)
            {
            p = v[j][i];
            v[j][i] = v[j][k];
            v[j][k] = p;
            }
         }
      }
}



/*---------------------------------------------------------------------------
   choldc() constructs the Cholesky decomposition of a given positive symmetric matrix

	input	
--------------------------------------------------------------------------*/
void choldc(double **amtx, int n, double *pvect)
{
   int i,j,k;
   double sum;
   
/*---
   for(i=1;i<=n;i++)
      {
      for(j=1;j<=n;j++)
         {
         printf("i=%d \t j=%d \t amtx=%f \n",i,j,amtx[i][j]);
         }
      } 
---*/ 
  
   for(i=1;i<=n;i++)
      {
      for(j=i;j<=n;j++)
         {
         sum=amtx[i][j];
         for(k=i-1;k>=1;k--)
            {
            sum -= amtx[i][k]*amtx[j][k];
            }
         if(i==j)
            {
            if(sum < 0.0)
               {
               printf("error in Cholesky decomposition \n");
               for(k=i;k>=1;k--) printf("i=%d \t k=%d \t amtx=%e \n",i,k,amtx[i][k]);
               getchar();
               exit(1);
               }
            pvect[i]=sqrt(sum);
            }
         else
            {
            amtx[j][i]=sum/pvect[i];
            }
         }
      }
}



/*---------------------------------------------------------------------------
   mult() gives the products of two matrices

	input	a	matrix a
		b	matrix b
		n	dimension of a and b
	output	c	matrix c = ab
---------------------------------------------------------------------------*/
void mult(double **a,double **b,int n,double **c)	
{
   int i,j,k;
   
   
   
   for(i=1;i<=n;i++)
      {
      for(j=1;j<=n;j++)
         {
         c[i][j]=0.0;
         for(k=1;k<=n;k++) c[i][j] += a[i][k]*b[k][j];
         }
      }
}



/*---------------------------------------------------------------------------
   solv1() gives the lower triangle of matrix C with CB=A in the following conditions:

 	A symmetric matrix
	B upper triangle matrix. Only the lower elements of C are computed.

	input	matrix a
		matrix b
		dimension of a and b
	output	matrix c = ab
--------------------------------------------------------------------------*/
void solv1(double **a,double **b,int n,double **c)	
{
   int i,j,k;
   double sum;
   
   for(i=1;i<=n;i++)
      {
      for(j=1;j<=i;j++)
         {
         sum = 0.0;
         for(k=1;k<j;k++) sum -= c[i][k]*b[k][j];
         c[i][j]=(a[j][i] + sum)/b[j][j];
         }
      }
}



/*---------------------------------------------------------------------------
   solv2() gives the lower triangle of C with AC=B in the following conditions:

	A lower triangle matrix
	C symmetric matrix
	B only the lower triangle of B is known.

	input	a	matrix A
		b	matrix B
		n	dimension of A and B
	output	c	matrix c
----------------------------------------------------------------------------*/
void solv2(double **a,double **b,int n,double **c)	
{
   int i,j,k;
   double sum;
   
   for(i=1;i<=n;i++)
      {
      for(j=1;j<=i;j++)
         {
         sum = 0.0;
         for(k=1;k<i;k++) sum -= a[i][k]*c[j][k];
         c[i][j] = c[j][i] = (b[i][j] + sum)/a[i][i];
         }
      }
}



/*----------------------------------------------------------------------------------------
   diagoz() solves the diagonalization problem:

	H.y = e.S.y  
 
	input:	nmax	matrix rank
		ss	overlapp matrix S
		hh	matrix to be diagonalized H
	output:	dvect	the eigenenegies vector sorted in increasing order
		y	the corresponding eigenvectors sorted accordingly
------------------------------------------------------------------------------------------*/
void diagoz(int nmax, double **ss, double **hh, double *d000, double **y000) 
{
   void choldc();
   void solv1();  			   
   void solv2();
   void tred2();
   void tqli();			
   void eigsrt();

   double **hmtx;
   double **smtx;

   double **cmtx;
   double *pvect;

   double **qmtx000;
   double *e000;
   
   double sum;
   
   int i,j,k;
   
/*   
printf("=========================DIAGONALISATION========================\n");
printf("nmax=%d \n",nmax); 
*/
   hmtx = dmatrix(1,nmax,1,nmax);
   smtx = dmatrix(1,nmax,1,nmax);
   cmtx = dmatrix(1,nmax,1,nmax);

   qmtx000 = dmatrix(1,nmax,1,nmax);
   e000 = dvector(1,nmax);
   
   pvect = dvector(1,nmax);
   
/* initialisation */
   for(i=1;i<=nmax;i++)
      {
      pvect[i]=0.0;
      for(j=1;j<=nmax;j++)
         {
         hmtx[i][j]=0.0;
         smtx[i][j]=0.0;
         cmtx[i][j]=0.0;
         y000[i][j]=0.0;
         }
      }
   
   for(i=1;i<=nmax;i++)					//fill-in hamiltonian h
      {
      for(j=1;j<=nmax;j++)
         {
         smtx[i][j] = ss[i][j];
         hmtx[i][j] = hh[i][j]; 
         }
      }

      
//---   printf("Cholesky decomposition \n");		// Cholesky decomposition of smtx
   choldc(smtx,nmax,pvect);
//---   printf("end of Cholesky decomposition \n");
   
   for(i=1;i<=nmax;i++)
      {
      smtx[i][i] = pvect[i];
      for(j=1;j<=nmax;j++) smtx[i][j]=smtx[j][i];
      }
      
   solv1(hmtx,smtx,nmax,y000);  			/* Y L^T = H */   
   solv2(smtx,y000,nmax,cmtx);				/* L C = Y */

   for(i=1;i<=nmax;i++)					/*fill-in qmtx000*/
      {
      for(j=1;j<=nmax;j++)
         {
         qmtx000[i][j] = cmtx[i][j];
         }
      }

   tred2(qmtx000,nmax,d000,e000);			/* tridiagonalisation de cmtx */
   tqli(d000,e000,nmax,qmtx000);			/* diagonalisation de cmtx*/

   for(i=1;i<=nmax;i++)					/* vecteurs propres de H = colonnes de la matrice y */
      {
      for(j=nmax;j>=1;j--)
         {
         sum = 0.0;
         y000[j][i] = 0.0;
         for(k=nmax;k>j;k--) sum -= smtx[j][k] * y000[k][i];   
         y000[j][i] = (qmtx000[j][i] + sum) / smtx[j][j];
         }
      } 

   eigsrt(d000,y000,nmax);
   
/*----
   for(i=1;i<=nmax;i++)
      {
      printf("i=%d \t d000=%e \n",i,d000[i]); 
      }

printf("---------fin resultat de la 1ere DIAGO apres classement--------\n");

getchar();
----*/

/*
for(i=1;i<=nmax;i++)
   {
   for(j=1;j<=nmax;j++)
      {
      if(fabs(y000[i][j]) > 1.0e-5) printf(" i=%d \t j=%d \t y000=%e \n",i,j,y000[i][j]);
      }
   getchar();
   } 
*/
        
/*----------------------------deallocation--------------------------*/   
   free_dmatrix(hmtx,1,nmax,1,nmax);
   free_dmatrix(smtx,1,nmax,1,nmax);
   free_dmatrix(qmtx000,1,nmax,1,nmax);
   free_dmatrix(cmtx,1,nmax,1,nmax);
   free_dvector(pvect,1,nmax);
   free_dvector(e000,1,nmax);
}




/*----------------------------------------------------------------------------------------
detdiagoz() computes the diagonalization of the hamiltonian of the determinants :   
 
 input :  pmax : size of the total hamiltonian 
          edet : matrix of the total hamiltonian < D1 | H | D2 >
          
 output:  ehh :  eigenvalues of the total hamiltonian
          yhh : eigenvectors of the total hamiltonian
------------------------------------------------------------------------------------------*/
void diagrs(int pmax, double **edet, double *ehh, double **yhh) 
{
   void tred2();  
   void tqli();	   			
   void eigsrt();

   double *evect;
   int i,j;
   
/*---   
   printf("=========================DIAGONALISATION========================\n");
   printf("Size of the total matrix pmax=%d \n",pmax); 
---*/
   evect = dvector(1,pmax);
   
   for(i=1;i<=pmax;i++)						
      {
      evect[i]=0.0;
      for(j=1;j<=pmax;j++)
         {
         yhh[i][j] = edet[i][j]; 
         }
      }
     
   tred2(yhh,pmax,ehh,evect); 
   tqli(ehh,evect,pmax,yhh);			
   eigsrt(ehh,yhh,pmax);

/*----------------------------deallocation--------------------------*/   
   free_dvector(evect,1,pmax);
}


void diag(int pmax, double **edet, double *ehh, double **yhh) 
{
   void tred2();  
   void tqli();	   			
   void eigsrt();

   double *evect;
   int i,j;
   
/*---   
   printf("=========================DIAGONALISATION========================\n");
   printf("Size of the total matrix pmax=%d \n",pmax); 
---*/
   evect = dvector(1,pmax);
   
   for(i=1;i<=pmax;i++)						
      {
      evect[i]=0.0;
      for(j=1;j<=pmax;j++)
         {
         yhh[i][j] = edet[i][j]; 
         }
      }
     
   tred2(yhh,pmax,ehh,evect); 
   tqli(ehh,evect,pmax,yhh);			
   eigsrt(ehh,yhh,pmax);

/*----------------------------deallocation--------------------------*/   
   free_dvector(evect,1,pmax);
}



/*----------------------------------------------------------------------------------------
   diagherm() computes the diagonalization of a complex hermitian matrix.
   After diagonalization, the yy vectors are sorted as follows:

		...	Y_2k-1	Y_2k	...
		---------------------------	
		...	-iV_k	 U_k	...
		...	  Uk	iV_k	...

 
	input :	nmax : matrix size 
		xxr : real part of the matrix
		xxi : imaginary part of the matrix
          
	output:	dx :  eigenvalues of the total hamiltonian
		qqr : real part of eigenvectors (U_k)
		qqi : imaginary part of eigenvectors (V_k)

   The frame transform is as follows:

           	 Xr = Qr.D.Qr^t + Qi.D.Qi^t
           	 Xi = Qi.D.Qr^t - Qr.D.Qi^t
   
------------------------------------------------------------------------------------------*/
void diagherm(int nmax, double **xxr, double **xxi, double **qqr, double **qqi, double *dx) 
{
   void tred2();  
   void tqli();	   			
   void eigsrt();
   double **yy;
   double *evect,*dvect;
   int nmax2;
   int i,j,i2,j2;

   nmax2 = 2*nmax;  
   evect = dvector(1,nmax2);
   dvect = dvector(1,nmax2);
   yy = dmatrix(1,nmax2,1,nmax2);

/*---   
   printf("=========================DIAGONALISATION========================\n");
   printf("Size of the total matrix nmax2=%d \n",nmax2); 
---*/
   
   for(i=1;i<=nmax;i++)				/*fill-in matrix yy*/
      {
      i2 = i + nmax;
      evect[i] = evect[i2] = dvect[i] = dvect[i2] = 0.0;
      for(j=1;j<=nmax;j++)
         {
         j2 = j + nmax;
         yy[i2][j2] = yy[i][j] = xxr[i][j]; 
         yy[i2][j] = xxi[i][j]; 
         yy[i][j2] = -xxi[i][j]; 
         }
      }
     
   tred2(yy,nmax2,dvect,evect); 
   tqli(dvect,evect,nmax2,yy);			
   eigsrt(dvect,yy,nmax2);

   for(i=1;i<=nmax;i++)
      {
      i2 = 2*i;
      dx[i] = dvect[i2];
      for(j=1;j<=nmax;j++)
         {
         j2 = j + nmax;
         qqr[j][i] = yy[j][i2]; 
         qqi[j][i] = yy[j2][i2]; 
         }
      }

/*---
   for(i=1;i<=nmax;i++)
      {
      i2 = i + nmax;
      for(j=1;j<=nmax;j++)
         {
         sum1 = sum2 = 0.0; 
         for(k=1;k<=nmax2;k++)
            {
            sum1 += dvect[k]*yy[i][k]*yy[j][k];
            sum2 += dvect[k]*yy[i2][k]*yy[j][k];
            } 
         printf("i=%d \t j=%d \t xxr=%e \t sum1=%e \t xxi=%e \t sum2=%e \n",i,j,xxr[i][j],sum1,xxi[i][j],sum2);
         }
getchar();
      }
---*/


/*---
   for(i=1;i<=nmax;i++)
      {
      for(j=1;j<=nmax;j++)
         {
         ttr[i][j] = tti[i][j] = 0.0; 
         for(k=1;k<=nmax;k++)
            {
            ttr[i][j] += dx[k]*(qqr[i][k]*qqr[j][k] + qqi[i][k]*qqi[j][k]);
            tti[i][j] += dx[k]*(qqi[i][k]*qqr[j][k] - qqr[i][k]*qqi[j][k]);
            } 
         printf("i=%d \t j=%d \t xxr=%e \t ttr=%e \t xxi=%e \t tti=%e \n",i,j,xxr[i][j],ttr[i][j],xxi[i][j],tti[i][j]);
         }
getchar();
      }
---*/


/*----------------------------deallocation--------------------------*/   
   free_dvector(evect,1,nmax2);
   free_dvector(dvect,1,nmax2);
   free_dmatrix(yy,1,nmax2,1,nmax2);
}
