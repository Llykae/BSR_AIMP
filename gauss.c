/*-------------------------------------------------------------------------
   File of mathematical constants and macros definitions
-------------------------------------------------------------------------*/

#include "gauss.h"

#define QMAX_GAUSS 48

double *egl;
double *wgl;

void init_gauss()
{
   switch(QMAX_GAUSS)
      {
      case 16:
	 egl = egl16;
	 wgl = wgl16;
	 break;

      case 24:
         egl = egl24;
         wgl = wgl24;
         break;

      case 48:
         egl = egl48;
         wgl = wgl48;
         break;

      default:
         printf("problem in init_gauss(): qmax_gauss=%d \n",QMAX_GAUSS);
         exit(1);
         break;
      }
}
