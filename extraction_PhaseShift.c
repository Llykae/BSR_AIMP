#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "nrutil.h"
#include "nrutil.c"

#define LMAX 2
#define file_length 42

void init(double *E, double *sigma, double **ratio);


void init(double *E, double *sigma, double **ratio)
{
    int k;
    double temp;

    FILE *fin;
    fin = fopen("input_extraction.data","r");
    for(k=1;k<=file_length;k++) fscanf(fin,"%lf %lf\n",&(E[k]),&(sigma[k]));
    fclose(fin);

    fin = fopen("./output/TCS/extraction/extraction_PS.out", "r");
    //for(k=1;k<=file_length;k++) fscanf(fin,"%lf %lf %lf %lf %lf\n",&temp,&(ratio[0][k]),&(ratio[1][k]),&(ratio[2][k]),&(ratio[3][k]));
    for(k=1;k<=file_length;k++) fscanf(fin,"%lf %lf %lf %lf \n",&temp,&(ratio[0][k]),&(ratio[1][k]),&(ratio[2][k]));
    fclose(fin);

    for(k=1;k<=file_length;k++) printf("%lf %lf %lf %lf %lf\n", E[k]*27.21, sigma[k], ratio[0][k], ratio[1][k], ratio[2][k]);
}


void main()
{
    int k, l;
    double sum;
    double *E, *sigma, **ratio;
    double **ps, *TCS;
    
    FILE *fout;
    fout = fopen("./output/TCS/extraction/ps_extracted.out", "w");
    
    E = dvector(1,file_length+1);
    sigma = dvector(1,file_length+1);
    ratio = dmatrix(0,LMAX,1,file_length+1);
    ps = dmatrix(0,LMAX,1,file_length+1);
    TCS = dvector(1,file_length+1);

    init(E,sigma,ratio);
    //getchar();

    for(k=1;k<=file_length;k++){
        sum=0.0;
        fprintf(fout,"%lf  \t",E[k]);
        for(l=0;l<=LMAX;l++){

            //ps[l][k] = asin(sqrt(E[k]*ratio[l][k]*(sigma[k]/(0.529*0.529))/(2*(2.*l+1.)*3.141592)));
            ps[l][k] = asin(sqrt(E[k]*ratio[l][k]*sigma[k]/(2*(2.*l+1.)*3.141592)));
            
            fprintf(fout," %lf  \t", ps[l][k]);

            sum+=(2.*l+1)*sin(ps[l][k])*sin(ps[l][k]);}
        fprintf(fout," %lf \t", 2*3.141592/E[k] * sum);
        fprintf(fout,"\n");}

    fclose(fout);   
}