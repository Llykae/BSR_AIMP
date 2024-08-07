/*--------------------------------------------------------------
   B-spline and related function evaluation:

		B. Gervais --- 2021
--------------------------------------------------------------*/
#ifndef BSP_UTILS_H
#define BSP_UTILS_H



/*--------------------------------------------------------------
   bsp_basis structures dexcribes the basis built on bsplines.
   Each basis function is made of a linear combination of n
   B-splines, defined by their index ib[1...n] and associated
   coefficients cb[1...n].
   The intervalle of definition of each basis function is given
   by [tmin,tmax]. The basis function vanish outside this
   intervalle
--------------------------------------------------------------*/
struct bsp_basis
{
   int ib;			//the unique index of the splines
//   int n;			//number of splines in the basis function
//   int *ib;			//the array of index
//   double *cb;			//the array of coefficients
   double tmin,tmax;		//basis function intervalle
};


double bspline(double *tn, int i, int k, double t);
double dbspline(double *tn, int i, int k, double t);
void bs_test(int nbs, int kb,  double *xi);
void basis_test(int nbs, int kb, double *xi, struct bsp_basis *bss);

#endif
