/*-------------------------------------------------------------------------------------------------
	This file contains the declaration of the function located gto.c
-------------------------------------------------------------------------------------------------*/
#ifndef PSEUD_EAU_H
#define PSEUD_EAU_H

double gtonorm(double a, int px, int py, int pz);
double fnorm(double a, int p1, int p2, int p3);
double gtoovrlp(struct orbital *bss, int i, int j, int k, int h);
double g1ovrlp(double mu, double nu, int m, int n, double xa, double xb);
void ovrlp_plnmlcoef(int m, int n, double xpa, double xpb, double *cp);
#endif
