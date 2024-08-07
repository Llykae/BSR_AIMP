#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*

void nerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   fprintf(stderr,"Numerical Recipes run-time error...\n") ;
   fprintf(stderr,"%s\n",error_text) ;
   fprintf(stderr,"...now exiting to system...\n") ;
   exit(1) ;
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
   float *v ;
   
   v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
   if (!v) nerror("allocation failure in vector()") ;
   return v-nl+NR_END ;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
   free((FREE_ARG) (v+nl-NR_END)) ;
}

double *dvector(long nl,long nh)
/*allocate a double vector with subscript range v[nl..nh]*/
{
   double *v ;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))) ;
   if (!v) printf ("allocation failure in dvector()") ;
   return v-nl+NR_END ;
}

void free_dvector(double *v, long nl, long nh)
/*free a double vector allocated with dvector()*/
{
   free((FREE_ARG) (v+nl-NR_END)) ;
}

int *ivector(long nl,long nh)
/*allocate a int vector with subscript range v[nl..nh]*/
{
   int *v ;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))) ;
   if (!v) printf ("allocation failure in ivector()") ;
   return v-nl+NR_END ;
}

void free_ivector(int *v, long nl, long nh)
/*free a int vector allocated with ivector()*/
{
   free((FREE_ARG) (v+nl-NR_END)) ;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
   long i, nrow=nrh-nrl+1, ncol=nch-ncl+1 ;
   float **m ;
   
   /* allocate pointers to rows */
   m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*))) ;

   if(!m)
   {
      printf("allocation failure 1 in matrix()") ;
      getchar() ;
   }
   m+=NR_END ;
   m-=nrl ;
   
   /* allocate rows and set pointers to them */
   m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float))) ;

   if (!m[nrl])
   {
      printf("allocation failure 2 in matrix()") ;
      getchar() ;
   }
   m[nrl]+=NR_END ;
   m[nrl]-=ncl ;
   
   for (i=nrl+1 ; i<=nrh ; i++) m[i]=m[i-1]+ncol ;
   
   /* return poiter to array of pointer to rows */
   return m ;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/*free a float matrix allocated with dmatrix()*/
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END)) ;
   free((FREE_ARG) (m+nrl-NR_END)) ;
}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
   long i, nrow=nrh-nrl+1, ncol=nch-ncl+1 ;
   double **m ;
   
   /* allocate pointers to rows */
   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*))) ;

   if(!m)
   {
      printf("allocation failure 1 in matrix()") ;
      getchar() ;
   }
   m+=NR_END ;
   m-=nrl ;
   
   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double))) ;

   if (!m[nrl])
   {
      printf("allocation failure 2 in matrix()") ;
      getchar() ;
   }
   m[nrl]+=NR_END ;
   m[nrl]-=ncl ;
   
   for (i=nrl+1 ; i<=nrh ; i++) m[i]=m[i-1]+ncol ;
   
   /* return poiter to array of pointer to rows */
   return m ;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/*free a double matrix allocated with dmatrix()*/
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END)) ;
   free((FREE_ARG) (m+nrl-NR_END)) ;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
   long i, nrow=nrh-nrl+1, ncol=nch-ncl+1 ;
   int **m ;
   
   /* allocate pointers to rows */
   m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*))) ;
   if (!m) 
   {
      printf("allocation failure 1 in matrix()") ;
      getchar() ;
   }
   m+=NR_END ;
   m-=nrl ;
   
   /* allocate rows and set pointers to them */
   m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int))) ;
   if (!m[nrl])
   {
      printf("allocation failure 2 in matrix()") ;
      getchar() ;
   }
   m[nrl]+=NR_END ;
   m[nrl]-=ncl ;
   
   for (i=nrl+1 ; i<=nrh ; i++) m[i]=m[i-1]+ncol ;
   
   /* return poiter to array of pointer to rows */
   return m ;
} 

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/*free a int matrix allocated with imatrix()*/
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END)) ;
   free((FREE_ARG) (m+nrl-NR_END)) ;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
   unsigned long *v;
   
   v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
   if (!v) nerror("allocation failure in lvector()");
   return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}

double ***f3tensor(long nrl, long nrh, long  ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
 long i,j,nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
 double ***t;
 
 /* allocate pointers to pointers to rows */
 t=(double ***) malloc((size_t)((nrow +NR_END)*sizeof(double**)));
 if (!t) nerror("allocation failor 1 in f3tensor()");
 t +=NR_END;
 t -=nrl;
 
 /* allocate pointers to rows and set pointers to them */
 t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
 if (!t[nrl]) nerror("allocation failure 2 in f3tensor()");
 t[nrl] +=NR_END;
 t[nrl] -= ncl;
 
 /*allocate rows and sey pointers to them */
 t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
 if (!t[nrl][ncl]) nerror("allocation failure 3 in f3tensor()"); 
 t[nrl][ncl] += NR_END;
 t[nrl][ncl] -= ndl;
 for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
 for(i=nrl+1;i<=nrh;i++) {
 	t[i]=t[i-1]+ncol;
 	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
 	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
 }
 
 return t;
}

void free_f3tensor(double ***t,long nrl, long nrh, long  ncl, long nch, long ndl, long ndh)
{ 
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}


