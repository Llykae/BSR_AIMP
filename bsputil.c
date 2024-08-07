/*--------------------------------------------------------------
   B-spline and related function evaluation:

	bsplines()

                        B. Gervais --- 2012
--------------------------------------------------------------*/

double bspline( double *tn, int i, int k, double t);
double dbspline( double *tn, int i, int k, double t);
void bs_test(int nb, int kb,  double *xi);

/*--------------------------------------------------------------
   bsplines evaluates the basis spline function B(i,k,x) at x
   for degree k and basis function i. recursive evaluation is
   done here.
--------------------------------------------------------------*/
double bspline( double *tn, int i, int k, double t)
{
   double b;
   double wi,wj;
   int j;


   j = i+1;
   if(k == 0)
      {
      if(t>=tn[i] && t <tn[j]) b = 1.0;
      else b = 0.0;
      }
   else
      {
      if(tn[i] == tn[i+k]) wi = 1.0;
      else
         {
         if(tn[i] < tn[i+k]) wi = (t-tn[i])/(tn[i+k] - tn[i]);
         else wi = 0.0;
         }
      if(tn[j] == tn[j+k]) wj =  1.0;
      else
         {
         if(tn[j] < tn[j+k]) wj = (t-tn[j])/(tn[j+k] - tn[j]);
         else wj = 0.0;
         }
      b = wi*bspline(tn,i,k-1,t) + (1.0 - wj)*bspline(tn,i+1,k-1,t);
      }

   return(b);
}


double dbspline( double *tn, int i, int k, double t)
{
   double db;
   //double wi,wj;
   int ik,j,jk,k1;

   j = i+1;
   ik = i + k;
   jk = j + k;
   k1 = k - 1;

   if(tn[i] == tn[jk])
      {
      db = 0.0;
      }
   else if((tn[i] < tn[ik]) && (tn[j] == tn[jk])) 
      {
      db = bspline(tn,i,k1,t)/(tn[ik] - tn[i]);
      db *= k;
      }
   else if((tn[i] == tn[ik]) && (tn[j] < tn[jk])) 
      {
      db = -bspline(tn,j,k1,t)/(tn[jk] - tn[j]);
      db *= k;
      }
   else 
      {
      db = bspline(tn,i,k1,t)/(tn[ik] - tn[i]) - bspline(tn,j,k1,t)/(tn[jk] - tn[j]);
      db *= k;
      }

   return(db);
}


void bs_test(int nb, int kb,  double *xi)
{
   FILE *fout;
   double x,bx[11],dbx[11];
   int j,m,n[11];

   fout = fopen("output/bsfunc.out","w");
   for(m=0;m<=10;m++) n[m] = 1 + m;

   for(j=0;j<=1000;j++)
      {
      x = xi[n[0]] + j*0.10/1000;
      fprintf(fout,"%e",x);
      for(m=0;m<=10;m++)
         {
         bx[m] = bspline(xi,n[m],kb,x);
         dbx[m] = dbspline(xi,n[m],kb,x);
         fprintf(fout,"\t%e\t%e",bx[m],dbx[m]);
         }
      fprintf(fout,"\n");
      }

   fclose(fout);
}



void outwave(int nb, int kb,  double *xi, double **corb)
{
   FILE *fout;
   double x,bkx,phimx;
   int j,k,m;//n[11];

   fout = fopen("output/eigenfunc.out","w");

   for(j=0;j<=1000;j++)
      {
      x = xi[1] + j*50.0/1000;
      fprintf(fout,"%e",x);
      for(m=1;m<=5;m++)
         {
         phimx = 0.0;
         for(k=1;k<=nb;k++)
            {
            bkx = bspline(xi,k,kb,x);
            phimx += corb[k][m]*bkx;
            }
         fprintf(fout,"\t%e",phimx);
         }
      fprintf(fout,"\n");
      }

   fclose(fout);
}

