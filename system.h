#ifndef __SYSTEM
#define __SYSTEM

    int ICONTINUUM = 1			//set continuum switch to ensure outer boundary condition
    int NOBLE_ATOM = 2            // 0 = He / 1 = Ne / 2 = Ar / 3 = Kr / 4 = Xe
    int LANG = 0					//angular momentum
    int LMAX = 3
    int Nk = 500
    int DKH = 1                   //SIGN = 0, Without DKH. SIGN = 1, With DKH.
    int TDSWC = 0                 //0 = TCS  -   1 = DCS

#endif