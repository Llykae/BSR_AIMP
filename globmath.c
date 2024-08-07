/*-------------------------------------------------------------------------
   File of mathematical constants and macros definitions
-------------------------------------------------------------------------*/
#define NLFAC 300
#define NBIN 150
#define NFAC 32
#define NBXMG 30
#define NBYMG 30

double pi=3.1415926535897932384;
double Pi=3.1415926535897932384;
double pi2=6.2831853071795864769;
double sqrtpi;
double sqpi;
double c=137.036;				//light velocity
double c2;					//light velocity squared


double lfac[NLFAC];				//log of the factorial
double dblfac[NFAC+1];				//double factorial : (n-1)!!
double fac[NFAC+1];				//first factorials
double bin[NBIN][NBIN];				//binomial coefficient
double xmg[NBXMG][NBXMG];			//xmg coefficient for use in hseries()
double ymg[NBYMG][NBYMG];			//ymg coefficient for use in hseries()
double zpl[NBYMG][NBYMG];			//zpl coefficient for use in uvseries()
double zpl1[NBYMG][NBYMG];			//zpl1 coefficient for use in uvseries()
double nrmzz[NFAC][NFAC];			//normalization coefficients for real spherical hrmonics

//-------------sum of reciprocal powers (zeta) from 1 to 42------------
double zeta[43]={0.00000000000000000000, 0.00000000000000000000,
                 1.64493406684822643647, 1.20205690315959428540,
                 1.08232323371113819152, 1.03692775514336992633,
                 1.01734306198444913971, 1.00834927738192282684,
                 1.00407735619794433938, 1.00200839282608221442,
                 1.00099457512781808534, 1.00049418860411946456,
                 1.00024608655330804830, 1.00012271334757848915,
                 1.00006124813505870483, 1.00003058823630702049,
                 1.00001528225940865187, 1.00000763719763789976,
                 1.00000381729326499984, 1.00000190821271655394,
                 1.00000095396203387280, 1.00000047693298678781,
                 1.00000023845050272773, 1.00000011921992596531,
                 1.00000005960818905126, 1.00000002980350351465,
                 1.00000001490155482837, 1.00000000745071178984,
                 1.00000000372533402479, 1.00000000186265972351,
                 1.00000000093132743242, 1.00000000046566290650,
                 1.00000000023283118337, 1.00000000011641550173,
                 1.00000000005820772088, 1.00000000002910385044,
                 1.00000000001455192189, 1.00000000000727595984,
                 1.00000000000363797955, 1.00000000000181898965,
                 1.00000000000090949478, 1.00000000000045474738,
                 1.00000000000022737368};


/*#define SIGN(a,b) ((b)>=0  ?  fabs(a) : -fabs(a))*/
#define MIN(a,b)  ((a)<(b)  ?  (a):(b))
#define MAX(a,b)  ((a)>(b)  ?  (a):(b))
#define SWAP(a,b) {double tmp; tmp;tmp=(a);(a)=(b);(b)=tmp;}  

