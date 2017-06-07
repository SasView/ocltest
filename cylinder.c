#line 1 "/home/pkienzle/src/sasmodels/sasmodels/kernel_header.c"
#ifdef __OPENCL_VERSION__
# define USE_OPENCL
#elif defined(_OPENMP)
# define USE_OPENMP
#endif

// If opencl is not available, then we are compiling a C function
// Note: if using a C++ compiler, then define kernel as extern "C"
#ifdef USE_OPENCL
   typedef int int32_t;
#  if defined(USE_SINCOS)
#    define SINCOS(angle,svar,cvar) svar=sincos(angle,&cvar)
#  else
#    define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
#  endif
   // Intel CPU on Mac gives strange values for erf(); on the verified
   // platforms (intel, nvidia, amd), the cephes erf() is significantly
   // faster than that available in the native OpenCL.
   #define NEED_ERF
   // OpenCL only has type generic math
   #define expf exp
   #ifndef NEED_ERF
   #  define erff erf
   #  define erfcf erfc
   #endif
#else // !USE_OPENCL
// Use SAS_DOUBLE to force the use of double even for float kernels
#  define SAS_DOUBLE dou ## ble
#  ifdef __cplusplus
      #include <cstdio>
      #include <cmath>
      using namespace std;
      #if defined(_MSC_VER)
         #include <limits>
         #include <float.h>
         #define kernel extern "C" __declspec( dllexport )
         inline double trunc(double x) { return x>=0?floor(x):-floor(-x); }
         inline double fmin(double x, double y) { return x>y ? y : x; }
         inline double fmax(double x, double y) { return x<y ? y : x; }
         #define isnan(x) _isnan(x)
         #define isinf(x) (!_finite(x))
         #define isfinite(x) _finite(x)
         #define NAN (std::numeric_limits<double>::quiet_NaN()) // non-signalling NaN
         #define INFINITY (std::numeric_limits<double>::infinity())
         #define NEED_ERF
         #define NEED_EXPM1
         #define NEED_TGAMMA
     #else
         #define kernel extern "C"
         #include <cstdint>
     #endif
     inline void SINCOS(double angle, double &svar, double &cvar) { svar=sin(angle); cvar=cos(angle); }
#  else // !__cplusplus
     #include <inttypes.h>  // C99 guarantees that int32_t types is here
     #include <stdio.h>
     #if defined(__TINYC__)
         typedef int int32_t;
         #include <math.h>
         // TODO: check isnan is correct
         inline double _isnan(double x) { return x != x; } // hope this doesn't optimize away!
         #undef isnan
         #define isnan(x) _isnan(x)
         // Defeat the double->float conversion since we don't have tgmath
         inline SAS_DOUBLE trunc(SAS_DOUBLE x) { return x>=0?floor(x):-floor(-x); }
         inline SAS_DOUBLE fmin(SAS_DOUBLE x, SAS_DOUBLE y) { return x>y ? y : x; }
         inline SAS_DOUBLE fmax(SAS_DOUBLE x, SAS_DOUBLE y) { return x<y ? y : x; }
         #define NEED_ERF
         #define NEED_EXPM1
         #define NEED_TGAMMA
         // expf missing from windows?
         #define expf exp
     #else
         #include <tgmath.h> // C99 type-generic math, so sin(float) => sinf
     #endif
     // MSVC doesn't support C99, so no need for dllexport on C99 branch
     #define kernel
     #define SINCOS(angle,svar,cvar) do {const double _t_=angle; svar=sin(_t_);cvar=cos(_t_);} while (0)
#  endif  // !__cplusplus
#  define global
#  define local
#  define constant const
// OpenCL powr(a,b) = C99 pow(a,b), b >= 0
// OpenCL pown(a,b) = C99 pow(a,b), b integer
#  define powr(a,b) pow(a,b)
#  define pown(a,b) pow(a,b)
#endif // !USE_OPENCL

#if defined(NEED_EXPM1)
   static SAS_DOUBLE expm1(SAS_DOUBLE x_in) {
      double x = (double)x_in;  // go back to float for single precision kernels
      // Adapted from the cephes math library.
      // Copyright 1984 - 1992 by Stephen L. Moshier
      if (x != x || x == 0.0) {
         return x; // NaN and +/- 0
      } else if (x < -0.5 || x > 0.5) {
         return exp(x) - 1.0;
      } else {
         const double xsq = x*x;
         const double p = (((
            +1.2617719307481059087798E-4)*xsq
            +3.0299440770744196129956E-2)*xsq
            +9.9999999999999999991025E-1);
         const double q = ((((
            +3.0019850513866445504159E-6)*xsq
            +2.5244834034968410419224E-3)*xsq
            +2.2726554820815502876593E-1)*xsq
            +2.0000000000000000000897E0);
         double r = x * p;
         r =  r / (q - r);
         return r+r;
       }
   }
#endif

// Standard mathematical constants:
//   M_E, M_LOG2E, M_LOG10E, M_LN2, M_LN10, M_PI, M_PI_2=pi/2, M_PI_4=pi/4,
//   M_1_PI=1/pi, M_2_PI=2/pi, M_2_SQRTPI=2/sqrt(pi), SQRT2, SQRT1_2=sqrt(1/2)
// OpenCL defines M_constant_F for float constants, and nothing if double
// is not enabled on the card, which is why these constants may be missing
#ifndef M_PI
#  define M_PI 3.141592653589793
#endif
#ifndef M_PI_2
#  define M_PI_2 1.570796326794897
#endif
#ifndef M_PI_4
#  define M_PI_4 0.7853981633974483
#endif
#ifndef M_E
#  define M_E 2.718281828459045091
#endif
#ifndef M_SQRT1_2
#  define M_SQRT1_2 0.70710678118654746
#endif

// Non-standard function library
// pi/180, used for converting between degrees and radians
// 4/3 pi for computing sphere volumes
// square and cube for computing squares and cubes
#ifndef M_PI_180
#  define M_PI_180 0.017453292519943295
#endif
#ifndef M_4PI_3
#  define M_4PI_3 4.18879020478639
#endif
inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }
inline double sas_sinx_x(double x) { return x==0 ? 1.0 : sin(x)/x; }

// To rotate from the canonical position to theta, phi, psi, first rotate by
// psi about the major axis, oriented along z, which is a rotation in the
// detector plane xy. Next rotate by theta about the y axis, aligning the major
// axis in the xz plane. Finally, rotate by phi in the detector plane xy.
// To compute the scattering, undo these rotations in reverse order:
//     rotate in xy by -phi, rotate in xz by -theta, rotate in xy by -psi
// The returned q is the length of the q vector and (xhat, yhat, zhat) is a unit
// vector in the q direction.
// To change between counterclockwise and clockwise rotation, change the
// sign of phi and psi.

#if 1
//think cos(theta) should be sin(theta) in new coords, RKH 11Jan2017
#define ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sn, cn) do { \
    SINCOS(phi*M_PI_180, sn, cn); \
    q = sqrt(qx*qx + qy*qy); \
    cn  = (q==0. ? 1.0 : (cn*qx + sn*qy)/q * sin(theta*M_PI_180));  \
    sn = sqrt(1 - cn*cn); \
    } while (0)
#else
// SasView 3.x definition of orientation
#define ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sn, cn) do { \
    SINCOS(theta*M_PI_180, sn, cn); \
    q = sqrt(qx*qx + qy*qy);\
    cn = (q==0. ? 1.0 : (cn*cos(phi*M_PI_180)*qx + sn*qy)/q); \
    sn = sqrt(1 - cn*cn); \
    } while (0)
#endif

#if 1
#define ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, xhat, yhat, zhat) do { \
    q = sqrt(qx*qx + qy*qy); \
    const double qxhat = qx/q; \
    const double qyhat = qy/q; \
    double sin_theta, cos_theta; \
    double sin_phi, cos_phi; \
    double sin_psi, cos_psi; \
    SINCOS(theta*M_PI_180, sin_theta, cos_theta); \
    SINCOS(phi*M_PI_180, sin_phi, cos_phi); \
    SINCOS(psi*M_PI_180, sin_psi, cos_psi); \
    xhat = qxhat*(-sin_phi*sin_psi + cos_theta*cos_phi*cos_psi) \
         + qyhat*( cos_phi*sin_psi + cos_theta*sin_phi*cos_psi); \
    yhat = qxhat*(-sin_phi*cos_psi - cos_theta*cos_phi*sin_psi) \
         + qyhat*( cos_phi*cos_psi - cos_theta*sin_phi*sin_psi); \
    zhat = qxhat*(-sin_theta*cos_phi) \
         + qyhat*(-sin_theta*sin_phi); \
    } while (0)
#else
// SasView 3.x definition of orientation
#define ORIENT_ASYMMETRIC(qx, qy, theta, phi, psi, q, cos_alpha, cos_mu, cos_nu) do { \
    q = sqrt(qx*qx + qy*qy); \
    const double qxhat = qx/q; \
    const double qyhat = qy/q; \
    double sin_theta, cos_theta; \
    double sin_phi, cos_phi; \
    double sin_psi, cos_psi; \
    SINCOS(theta*M_PI_180, sin_theta, cos_theta); \
    SINCOS(phi*M_PI_180, sin_phi, cos_phi); \
    SINCOS(psi*M_PI_180, sin_psi, cos_psi); \
    cos_alpha = cos_theta*cos_phi*qxhat + sin_theta*qyhat; \
    cos_mu = (-sin_theta*cos_psi*cos_phi - sin_psi*sin_phi)*qxhat + cos_theta*cos_psi*qyhat; \
    cos_nu = (-cos_phi*sin_psi*sin_theta + sin_phi*cos_psi)*qxhat + sin_psi*cos_theta*qyhat; \
    } while (0)
#endif

#line 1 "/home/pkienzle/src/sasmodels/sasmodels/models/lib/polevl.c"
/*							polevl.c
 *							p1evl.c
 *
 *	Evaluate polynomial
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N+1], polevl[];
 *
 * y = polevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates polynomial of degree N:
 *
 *                     2          N
 * y  =  C  + C x + C x  +...+ C x
 *        0    1     2          N
 *
 * Coefficients are stored in reverse order:
 *
 * coef[0] = C  , ..., coef[N] = C  .
 *            N                   0
 *
 * The function p1evl() assumes that C_N = 1.0 and is
 * omitted from the array.  Its calling arguments are
 * otherwise the same as polevl().
 *
 *
 * SPEED:
 *
 * In the interest of speed, there are no checks for out
 * of bounds arithmetic.  This routine is used by most of
 * the functions in the library.  Depending on available
 * equipment features, the user may wish to rewrite the
 * program in microcode or assembly language.
 *
 */


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

double polevl( double x, constant double *coef, int N );
double polevl( double x, constant double *coef, int N )
{

    int i = 0;
    double ans = coef[i];

    while (i < N) {
        i++;
        ans = ans * x + coef[i];
    }

    return ans;
}

/*							p1evl()	*/
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl( double x, constant double *coef, int N );
double p1evl( double x, constant double *coef, int N )
{
    int i=0;
    double ans = x+coef[i];

    while (i < N-1) {
        i++;
        ans = ans*x + coef[i];
    }

    return ans;
}
#line 1 "/home/pkienzle/src/sasmodels/sasmodels/models/lib/sas_J1.c"
/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

#if FLOAT_SIZE>4
//Cephes double pression function
double cephes_j1(double x);

constant double RPJ1[8] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
    0.0,
    0.0,
    0.0,
    0.0 };

constant double RQJ1[8] = {
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18
    };

constant double PPJ1[8] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
    0.0} ;


constant double PQJ1[8] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
    0.0 };

constant double QPJ1[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1 };

constant double QQJ1[8] = {
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
    0.0 };

double cephes_j1(double x)
{

    double w, z, p, q, abs_x, sign_x;

    const double Z1 = 1.46819706421238932572E1;
    const double Z2 = 4.92184563216946036703E1;

    // 2017-05-18 PAK - mathematica and mpmath use J1(-x) = -J1(x)
    if (x < 0) {
        abs_x = -x;
        sign_x = -1.0;
    } else {
        abs_x = x;
        sign_x = 1.0;
    }

    if( abs_x <= 5.0 ) {
        z = abs_x * abs_x;
        w = polevl( z, RPJ1, 3 ) / p1evl( z, RQJ1, 8 );
        w = w * abs_x * (z - Z1) * (z - Z2);
        return( sign_x * w );
    }

    w = 5.0/abs_x;
    z = w * w;
    p = polevl( z, PPJ1, 6)/polevl( z, PQJ1, 6 );
    q = polevl( z, QPJ1, 7)/p1evl( z, QQJ1, 7 );

    // 2017-05-19 PAK improve accuracy using trig identies
    // original:
    //    const double THPIO4 =  2.35619449019234492885;
    //    const double SQ2OPI = 0.79788456080286535588;
    //    double sin_xn, cos_xn;
    //    SINCOS(abs_x - THPIO4, sin_xn, cos_xn);
    //    p = p * cos_xn - w * q * sin_xn;
    //    return( sign_x * p * SQ2OPI / sqrt(abs_x) );
    // expanding p*cos(a - 3 pi/4) - wq sin(a - 3 pi/4)
    //    [ p(sin(a) - cos(a)) + wq(sin(a) + cos(a)) / sqrt(2)
    // note that sqrt(1/2) * sqrt(2/pi) = sqrt(1/pi)
    const double SQRT1_PI = 0.56418958354775628;
    double sin_x, cos_x;
    SINCOS(abs_x, sin_x, cos_x);
    p = p*(sin_x - cos_x) + w*q*(sin_x + cos_x);
    return( sign_x * p * SQRT1_PI / sqrt(abs_x) );
}

#else
//Single precission version of cephes
float cephes_j1f(float x);

constant float JPJ1[8] = {
    -4.878788132172128E-009,
    6.009061827883699E-007,
    -4.541343896997497E-005,
    1.937383947804541E-003,
    -3.405537384615824E-002,
    0.0,
    0.0,
    0.0
    };

constant float MO1J1[8] = {
    6.913942741265801E-002,
    -2.284801500053359E-001,
    3.138238455499697E-001,
    -2.102302420403875E-001,
    5.435364690523026E-003,
    1.493389585089498E-001,
    4.976029650847191E-006,
    7.978845453073848E-001
    };

constant float PH1J1[8] = {
    -4.497014141919556E+001,
    5.073465654089319E+001,
    -2.485774108720340E+001,
    7.222973196770240E+000,
    -1.544842782180211E+000,
    3.503787691653334E-001,
    -1.637986776941202E-001,
    3.749989509080821E-001
    };

float cephes_j1f(float xx)
{

    float x, w, z, p, q, xn;

    const float Z1 = 1.46819706421238932572E1;


    // 2017-05-18 PAK - mathematica and mpmath use J1(-x) = -J1(x)
    x = xx;
    if( x < 0 )
        x = -xx;

    if( x <= 2.0 ) {
        z = x * x;
        p = (z-Z1) * x * polevl( z, JPJ1, 4 );
        return( xx < 0. ? -p : p );
    }

    q = 1.0/x;
    w = sqrt(q);

    p = w * polevl( q, MO1J1, 7);
    w = q*q;
    // 2017-05-19 PAK improve accuracy using trig identies
    // original:
    //    const float THPIO4F =  2.35619449019234492885;    /* 3*pi/4 */
    //    xn = q * polevl( w, PH1J1, 7) - THPIO4F;
    //    p = p * cos(xn + x);
    //    return( xx < 0. ? -p : p );
    // expanding cos(a + b - 3 pi/4) is
    //    [sin(a)sin(b) + sin(a)cos(b) + cos(a)sin(b)-cos(a)cos(b)] / sqrt(2)
    xn = q * polevl( w, PH1J1, 7);
    float cos_xn, sin_xn;
    float cos_x, sin_x;
    SINCOS(xn, sin_xn, cos_xn);  // about xn and 1
    SINCOS(x, sin_x, cos_x);
    p *= M_SQRT1_2*(sin_xn*(sin_x+cos_x) + cos_xn*(sin_x-cos_x));

    return( xx < 0. ? -p : p );
}
#endif

#if FLOAT_SIZE>4
#define sas_J1 cephes_j1
#else
#define sas_J1 cephes_j1f
#endif

//Finally J1c function that equals 2*J1(x)/x
double sas_2J1x_x(double x);
double sas_2J1x_x(double x)
{
    return (x != 0.0 ) ? 2.0*sas_J1(x)/x : 1.0;
}
#line 1 "/home/pkienzle/src/sasmodels/sasmodels/models/lib/gauss76.c"
/*
 *  GaussWeights.c
 *  SANSAnalysis
 *
 *  Created by Andrew Jackson on 4/23/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#define N_POINTS_76 76

// Gaussians
constant double Gauss76Wt[N_POINTS_76]={
	.00126779163408536,		//0
	.00294910295364247,
	.00462793522803742,
	.00629918049732845,
	.00795984747723973,
	.00960710541471375,
	.0112381685696677,
	.0128502838475101,
	.0144407317482767,
	.0160068299122486,
	.0175459372914742,		//10
	.0190554584671906,
	.020532847967908,
	.0219756145344162,
	.0233813253070112,
	.0247476099206597,
	.026072164497986,
	.0273527555318275,
	.028587223650054,
	.029773487255905,
	.0309095460374916,		//20
	.0319934843404216,
	.0330234743977917,
	.0339977794120564,
	.0349147564835508,
	.0357728593807139,
	.0365706411473296,
	.0373067565423816,
	.0379799643084053,
	.0385891292645067,
	.0391332242205184,		//30
	.0396113317090621,
	.0400226455325968,
	.040366472122844,
	.0406422317102947,
	.0408494593018285,
	.040987805464794,
	.0410570369162294,
	.0410570369162294,
	.040987805464794,
	.0408494593018285,		//40
	.0406422317102947,
	.040366472122844,
	.0400226455325968,
	.0396113317090621,
	.0391332242205184,
	.0385891292645067,
	.0379799643084053,
	.0373067565423816,
	.0365706411473296,
	.0357728593807139,		//50
	.0349147564835508,
	.0339977794120564,
	.0330234743977917,
	.0319934843404216,
	.0309095460374916,
	.029773487255905,
	.028587223650054,
	.0273527555318275,
	.026072164497986,
	.0247476099206597,		//60
	.0233813253070112,
	.0219756145344162,
	.020532847967908,
	.0190554584671906,
	.0175459372914742,
	.0160068299122486,
	.0144407317482767,
	.0128502838475101,
	.0112381685696677,
	.00960710541471375,		//70
	.00795984747723973,
	.00629918049732845,
	.00462793522803742,
	.00294910295364247,
	.00126779163408536		//75 (indexed from 0)
};

constant double Gauss76Z[N_POINTS_76]={
	-.999505948362153,		//0
	-.997397786355355,
	-.993608772723527,
	-.988144453359837,
	-.981013938975656,
	-.972229228520377,
	-.961805126758768,
	-.949759207710896,
	-.936111781934811,
	-.92088586125215,
	-.904107119545567,		//10
	-.885803849292083,
	-.866006913771982,
	-.844749694983342,
	-.822068037328975,
	-.7980001871612,
	-.77258672828181,
	-.74587051350361,
	-.717896592387704,
	-.688712135277641,
	-.658366353758143,		//20
	-.626910417672267,
	-.594397368836793,
	-.560882031601237,
	-.526420920401243,
	-.491072144462194,
	-.454895309813726,
	-.417951418780327,
	-.380302767117504,
	-.342012838966962,
	-.303146199807908,		//30
	-.263768387584994,
	-.223945802196474,
	-.183745593528914,
	-.143235548227268,
	-.102483975391227,
	-.0615595913906112,
	-.0205314039939986,
	.0205314039939986,
	.0615595913906112,
	.102483975391227,			//40
	.143235548227268,
	.183745593528914,
	.223945802196474,
	.263768387584994,
	.303146199807908,
	.342012838966962,
	.380302767117504,
	.417951418780327,
	.454895309813726,
	.491072144462194,		//50
	.526420920401243,
	.560882031601237,
	.594397368836793,
	.626910417672267,
	.658366353758143,
	.688712135277641,
	.717896592387704,
	.74587051350361,
	.77258672828181,
	.7980001871612,	//60
	.822068037328975,
	.844749694983342,
	.866006913771982,
	.885803849292083,
	.904107119545567,
	.92088586125215,
	.936111781934811,
	.949759207710896,
	.961805126758768,
	.972229228520377,		//70
	.981013938975656,
	.988144453359837,
	.993608772723527,
	.997397786355355,
	.999505948362153		//75
};

#line 1 "/home/pkienzle/src/sasmodels/sasmodels/models/cylinder.c"
double form_volume(double radius, double length);
double fq(double q, double sn, double cn,double radius, double length);
double orient_avg_1D(double q, double radius, double length);
double Iq(double q, double sld, double solvent_sld, double radius, double length);
double Iqxy(double qx, double qy, double sld, double solvent_sld,
    double radius, double length, double theta, double phi);

#define INVALID(v) (v.radius<0 || v.length<0)

double form_volume(double radius, double length)
{
    return M_PI*radius*radius*length;
}

double fq(double q, double sn, double cn, double radius, double length)
{
    // precompute qr and qh to save time in the loop
    const double qr = q*radius;
    const double qh = q*0.5*length; 
    return sas_2J1x_x(qr*sn) * sas_sinx_x(qh*cn);
}

double orient_avg_1D(double q, double radius, double length)
{
    // translate a point in [-1,1] to a point in [0, pi/2]
    const double zm = M_PI_4;
    const double zb = M_PI_4; 

    double total = 0.0;
    for (int i=0; i<76 ;i++) {
        const double alpha = Gauss76Z[i]*zm + zb;
        double sn, cn; // slots to hold sincos function output
        // alpha(theta,phi) the projection of the cylinder on the detector plane
        SINCOS(alpha, sn, cn);
        total += Gauss76Wt[i] * square( fq(q, sn, cn, radius, length) ) * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    return total*zm;
}

double Iq(double q,
    double sld,
    double solvent_sld,
    double radius,
    double length)
{
    const double s = (sld - solvent_sld) * form_volume(radius, length);
    return 1.0e-4 * s * s * orient_avg_1D(q, radius, length);
}


double Iqxy(double qx, double qy,
    double sld,
    double solvent_sld,
    double radius,
    double length,
    double theta,
    double phi)
{
    double q, sin_alpha, cos_alpha;
    ORIENT_SYMMETRIC(qx, qy, theta, phi, q, sin_alpha, cos_alpha);
    //printf("sld: %g solvent_sld: %g radius: %g length: %g theta: %g phi: %g\n", sld, solvent_sld, radius, length, theta, phi);
    //printf("qx: %g qy: %g theta:%g phi: %g sn: %g cn: %g\n", qx, qy, theta, phi, sin_alpha, cos_alpha);
    const double s = (sld-solvent_sld) * form_volume(radius, length);
    const double form = fq(q, sin_alpha, cos_alpha, radius, length);
    return 1.0e-4 * square(s * form);
}

#line 542 "sasmodels/generate.py"
#define PARAMETER_TABLE \
double sld;\
double sld_solvent;\
double radius;\
double length;\
double theta;\
double phi;
#define CALL_VOLUME(_v) form_volume(_v.radius,_v.length)
#define MAX_PD 4
#define NUM_PARS 6
#define NUM_VALUES 17
#define NUM_MAGNETIC 2
#define MAGNETIC_PARS 0,1
#define MAGNETIC_PAR1 0
#define MAGNETIC_PAR2 1

#if 0 // cylinder Iq
#define KERNEL_NAME cylinder_Iq
#define CALL_IQ(_q,_i,_v) Iq(_q[_i],_v.sld,_v.sld_solvent,_v.radius,_v.length)
#line 1 "/home/pkienzle/src/sasmodels/sasmodels/kernel_iq.cl Iq"

/*
    ##########################################################
    #                                                        #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                              !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNELPY.PY  !!   #
    #   !!                                              !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                        #
    ##########################################################
*/

#ifndef _PAR_BLOCK_ // protected block so we can include this code twice.
#define _PAR_BLOCK_

typedef struct {
#if MAX_PD > 0
    int32_t pd_par[MAX_PD];     // id of the nth polydispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the value & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
#endif // MAX_PD > 0
    int32_t num_eval;           // total number of voxels in hypercube
    int32_t num_weights;        // total length of the weights vector
    int32_t num_active;         // number of non-trivial pd loops
    int32_t theta_par;          // id of spherical correction variable
} ProblemDetails;

// Intel HD 4000 needs private arrays to be a multiple of 4 long
typedef struct {
    PARAMETER_TABLE
} ParameterTable;
typedef union {
    ParameterTable table;
    double vector[4*((NUM_PARS+3)/4)];
} ParameterBlock;
#endif // _PAR_BLOCK_


#if defined(MAGNETIC) && NUM_MAGNETIC>0

// Return value restricted between low and high
static double clip(double value, double low, double high)
{
  return (value < low ? low : (value > high ? high : value));
}

// Compute spin cross sections given in_spin and out_spin
// To convert spin cross sections to sld b:
//     uu * (sld - m_sigma_x);
//     dd * (sld + m_sigma_x);
//     ud * (m_sigma_y + 1j*m_sigma_z);
//     du * (m_sigma_y - 1j*m_sigma_z);
static void set_spins(double in_spin, double out_spin, double spins[4])
{
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  spins[0] = sqrt(sqrt((1.0-in_spin) * (1.0-out_spin))); // dd
  spins[1] = sqrt(sqrt((1.0-in_spin) * out_spin));       // du
  spins[2] = sqrt(sqrt(in_spin * (1.0-out_spin)));       // ud
  spins[3] = sqrt(sqrt(in_spin * out_spin));             // uu
}

static double mag_sld(double qx, double qy, double p,
                       double mx, double my, double sld)
{
    const double perp = qy*mx - qx*my;
    return sld + perp*p;
}

#endif // MAGNETIC

kernel
void KERNEL_NAME(
    int32_t nq,                 // number of q values
    const int32_t pd_start,     // where we are in the polydispersity loop
    const int32_t pd_stop,      // where we are stopping in the polydispersity loop
    global const ProblemDetails *details,
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+1 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{

  // who we are and what element we are working with
  const int q_index = get_global_id(0);
  if (q_index >= nq) return;

  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  ParameterBlock local_values;

#if defined(MAGNETIC) && NUM_MAGNETIC>0
  // Location of the sld parameters in the parameter vector.
  // These parameters are updated with the effective sld due to magnetism.
  #if NUM_MAGNETIC > 3
  const int32_t slds[] = { MAGNETIC_PARS };
  #endif

  // TODO: could precompute these outside of the kernel.
  // Interpret polarization cross section.
  //     up_frac_i = values[NUM_PARS+2];
  //     up_frac_f = values[NUM_PARS+3];
  //     up_angle = values[NUM_PARS+4];
  double spins[4];
  double cos_mspin, sin_mspin;
  set_spins(values[NUM_PARS+2], values[NUM_PARS+3], spins);
  SINCOS(-values[NUM_PARS+4]*M_PI_180, sin_mspin, cos_mspin);
#endif // MAGNETIC

  // Fill in the initial variables
  //   values[0] is scale
  //   values[1] is background
  for (int i=0; i < NUM_PARS; i++) {
    local_values.vector[i] = values[2+i];
//if (q_index==0) printf("p%d = %g\n",i, local_values.vector[i]);
  }
//if (q_index==0) printf("NUM_VALUES:%d  NUM_PARS:%d  MAX_PD:%d\n", NUM_VALUES, NUM_PARS, MAX_PD);
//if (q_index==0) printf("start:%d stop:%d\n", pd_start, pd_stop);

  double pd_norm = (pd_start == 0 ? 0.0 : result[nq]);
  double this_result = (pd_start == 0 ? 0.0 : result[q_index]);
//if (q_index==0) printf("start %d %g %g\n", pd_start, pd_norm, this_result);

#if MAX_PD>0
  global const double *pd_value = values + NUM_VALUES;
  global const double *pd_weight = pd_value + details->num_weights;
#endif

  // Jump into the middle of the polydispersity loop
#if MAX_PD>4
  int n4=details->pd_length[4];
  int i4=(pd_start/details->pd_stride[4])%n4;
  const int p4=details->pd_par[4];
  global const double *v4 = pd_value + details->pd_offset[4];
  global const double *w4 = pd_weight + details->pd_offset[4];
#endif
#if MAX_PD>3
  int n3=details->pd_length[3];
  int i3=(pd_start/details->pd_stride[3])%n3;
  const int p3=details->pd_par[3];
  global const double *v3 = pd_value + details->pd_offset[3];
  global const double *w3 = pd_weight + details->pd_offset[3];
//if (q_index==0) printf("offset %d: %d %d\n", 3, details->pd_offset[3], NUM_VALUES);
#endif
#if MAX_PD>2
  int n2=details->pd_length[2];
  int i2=(pd_start/details->pd_stride[2])%n2;
  const int p2=details->pd_par[2];
  global const double *v2 = pd_value + details->pd_offset[2];
  global const double *w2 = pd_weight + details->pd_offset[2];
#endif
#if MAX_PD>1
  int n1=details->pd_length[1];
  int i1=(pd_start/details->pd_stride[1])%n1;
  const int p1=details->pd_par[1];
  global const double *v1 = pd_value + details->pd_offset[1];
  global const double *w1 = pd_weight + details->pd_offset[1];
#endif
#if MAX_PD>0
  int n0=details->pd_length[0];
  int i0=(pd_start/details->pd_stride[0])%n0;
  const int p0=details->pd_par[0];
  global const double *v0 = pd_value + details->pd_offset[0];
  global const double *w0 = pd_weight + details->pd_offset[0];
#endif


#if MAX_PD>0
  const int theta_par = details->theta_par;
  const bool fast_theta = (theta_par == p0);
  const bool slow_theta = (theta_par >= 0 && !fast_theta);
  double spherical_correction = 1.0;
#else
  // Note: if not polydisperse the weights cancel and we don't need the
  // spherical correction.
  const double spherical_correction = 1.0;
#endif

  int step = pd_start;


#if MAX_PD>4
  const double weight5 = 1.0;
  while (i4 < n4) {
    local_values.vector[p4] = v4[i4];
    double weight4 = w4[i4] * weight5;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 4, p4, i4, n4, local_values.vector[p4], weight4);
#elif MAX_PD>3
    const double weight4 = 1.0;
#endif
#if MAX_PD>3
  while (i3 < n3) {
    local_values.vector[p3] = v3[i3];
    double weight3 = w3[i3] * weight4;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 3, p3, i3, n3, local_values.vector[p3], weight3);
#elif MAX_PD>2
    const double weight3 = 1.0;
#endif
#if MAX_PD>2
  while (i2 < n2) {
    local_values.vector[p2] = v2[i2];
    double weight2 = w2[i2] * weight3;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 2, p2, i2, n2, local_values.vector[p2], weight2);
#elif MAX_PD>1
    const double weight2 = 1.0;
#endif
#if MAX_PD>1
  while (i1 < n1) {
    local_values.vector[p1] = v1[i1];
    double weight1 = w1[i1] * weight2;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 1, p1, i1, n1, local_values.vector[p1], weight1);
#elif MAX_PD>0
    const double weight1 = 1.0;
#endif
#if MAX_PD>0
  if (slow_theta) { // Theta is not in inner loop
    spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[theta_par])), 1.e-6);
  }
  while(i0 < n0) {
    local_values.vector[p0] = v0[i0];
    double weight0 = w0[i0] * weight1;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 0, p0, i0, n0, local_values.vector[p0], weight0);
    if (fast_theta) { // Theta is in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[p0])), 1.e-6);
    }
#else
    const double weight0 = 1.0;
#endif

//if (q_index == 0) {printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NUM_PARS; i++) printf("p%d=%g ",i, local_values.vector[i]); printf("\n"); }
//if (q_index == 0) printf("sphcor: %g\n", spherical_correction);

    #ifdef INVALID
    if (!INVALID(local_values.table))
    #endif
    {
      // Accumulate I(q)
      // Note: weight==0 must always be excluded
      if (weight0 > cutoff) {
        // spherical correction is set at a minimum of 1e-6, otherwise there
        // would be problems looking at models with theta=90.
        const double weight = weight0 * spherical_correction;
        pd_norm += weight * CALL_VOLUME(local_values.table);

#if defined(MAGNETIC) && NUM_MAGNETIC > 0
        const double qx = q[2*q_index];
        const double qy = q[2*q_index+1];
        const double qsq = qx*qx + qy*qy;

        // Constant across orientation, polydispersity for given qx, qy
        double scattering = 0.0;
        // TODO: what is the magnetic scattering at q=0
        if (qsq > 1.e-16) {
          double p[4];  // dd, du, ud, uu
          p[0] = (qy*cos_mspin + qx*sin_mspin)/qsq;
          p[3] = -p[0];
          p[1] = p[2] = (qy*sin_mspin - qx*cos_mspin)/qsq;

          for (int index=0; index<4; index++) {
            const double xs = spins[index];
            if (xs > 1.e-8) {
              const int spin_flip = (index==1) || (index==2);
              const double pk = p[index];
              for (int axis=0; axis<=spin_flip; axis++) {
                #define M1 NUM_PARS+5
                #define M2 NUM_PARS+8
                #define M3 NUM_PARS+13
                #define SLD(_M_offset, _sld_offset) \
                    local_values.vector[_sld_offset] = xs * (axis \
                    ? (index==1 ? -values[_M_offset+2] : values[_M_offset+2]) \
                    : mag_sld(qx, qy, pk, values[_M_offset], values[_M_offset+1], \
                              (spin_flip ? 0.0 : values[_sld_offset+2])))
                #if NUM_MAGNETIC==1
                    SLD(M1, MAGNETIC_PAR1);
                #elif NUM_MAGNETIC==2
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                #elif NUM_MAGNETIC==3
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                    SLD(M3, MAGNETIC_PAR3);
                #else
                for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                    SLD(M1+3*sk, slds[sk]);
                }
                #endif
                scattering += CALL_IQ(q, q_index, local_values.table);
              }
            }
          }
        }
#else  // !MAGNETIC
        const double scattering = CALL_IQ(q, q_index, local_values.table);
#endif // !MAGNETIC
        this_result += weight * scattering;
      }
    }
    ++step;
#if MAX_PD>0
    if (step >= pd_stop) break;
    ++i0;
  }
  i0 = 0;
#endif
#if MAX_PD>1
    if (step >= pd_stop) break;
    ++i1;
  }
  i1 = 0;
#endif
#if MAX_PD>2
    if (step >= pd_stop) break;
    ++i2;
  }
  i2 = 0;
#endif
#if MAX_PD>3
    if (step >= pd_stop) break;
    ++i3;
  }
  i3 = 0;
#endif
#if MAX_PD>4
    if (step >= pd_stop) break;
    ++i4;
  }
  i4 = 0;
#endif

//if (q_index==0) printf("res: %g/%g\n", this_result, pd_norm);
  // Remember the current result and the updated norm.
  result[q_index] = this_result;
  if (q_index == 0) result[nq] = pd_norm;
}

#undef CALL_IQ
#undef KERNEL_NAME
#endif // cylinder_Iq

#if 1 // cylinder_Iqxy
#define KERNEL_NAME cylinder_Iqxy
#define CALL_IQ(_q,_i,_v) Iqxy(_q[2*_i],_q[2*_i+1],_v.sld,_v.sld_solvent,_v.radius,_v.length,_v.theta,_v.phi)
#line 1 "/home/pkienzle/src/sasmodels/sasmodels/kernel_iq.cl Iqxy"

/*
    ##########################################################
    #                                                        #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                              !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNELPY.PY  !!   #
    #   !!                                              !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                        #
    ##########################################################
*/

#ifndef _PAR_BLOCK_ // protected block so we can include this code twice.
#define _PAR_BLOCK_

typedef struct {
#if MAX_PD > 0
    int32_t pd_par[MAX_PD];     // id of the nth polydispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the value & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
#endif // MAX_PD > 0
    int32_t num_eval;           // total number of voxels in hypercube
    int32_t num_weights;        // total length of the weights vector
    int32_t num_active;         // number of non-trivial pd loops
    int32_t theta_par;          // id of spherical correction variable
} ProblemDetails;

// Intel HD 4000 needs private arrays to be a multiple of 4 long
typedef struct {
    PARAMETER_TABLE
} ParameterTable;
typedef union {
    ParameterTable table;
    double vector[4*((NUM_PARS+3)/4)];
} ParameterBlock;
#endif // _PAR_BLOCK_


#if defined(MAGNETIC) && NUM_MAGNETIC>0

// Return value restricted between low and high
static double clip(double value, double low, double high)
{
  return (value < low ? low : (value > high ? high : value));
}

// Compute spin cross sections given in_spin and out_spin
// To convert spin cross sections to sld b:
//     uu * (sld - m_sigma_x);
//     dd * (sld + m_sigma_x);
//     ud * (m_sigma_y + 1j*m_sigma_z);
//     du * (m_sigma_y - 1j*m_sigma_z);
static void set_spins(double in_spin, double out_spin, double spins[4])
{
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  spins[0] = sqrt(sqrt((1.0-in_spin) * (1.0-out_spin))); // dd
  spins[1] = sqrt(sqrt((1.0-in_spin) * out_spin));       // du
  spins[2] = sqrt(sqrt(in_spin * (1.0-out_spin)));       // ud
  spins[3] = sqrt(sqrt(in_spin * out_spin));             // uu
}

static double mag_sld(double qx, double qy, double p,
                       double mx, double my, double sld)
{
    const double perp = qy*mx - qx*my;
    return sld + perp*p;
}

#endif // MAGNETIC

kernel
void KERNEL_NAME(
    int32_t nq,                 // number of q values
    const int32_t pd_start,     // where we are in the polydispersity loop
    const int32_t pd_stop,      // where we are stopping in the polydispersity loop
    global const ProblemDetails *details,
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+1 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{

  // who we are and what element we are working with
  const int q_index = get_global_id(0);
  if (q_index >= nq) return;

  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  ParameterBlock local_values;

#if defined(MAGNETIC) && NUM_MAGNETIC>0
  // Location of the sld parameters in the parameter vector.
  // These parameters are updated with the effective sld due to magnetism.
  #if NUM_MAGNETIC > 3
  const int32_t slds[] = { MAGNETIC_PARS };
  #endif

  // TODO: could precompute these outside of the kernel.
  // Interpret polarization cross section.
  //     up_frac_i = values[NUM_PARS+2];
  //     up_frac_f = values[NUM_PARS+3];
  //     up_angle = values[NUM_PARS+4];
  double spins[4];
  double cos_mspin, sin_mspin;
  set_spins(values[NUM_PARS+2], values[NUM_PARS+3], spins);
  SINCOS(-values[NUM_PARS+4]*M_PI_180, sin_mspin, cos_mspin);
#endif // MAGNETIC

  // Fill in the initial variables
  //   values[0] is scale
  //   values[1] is background
  for (int i=0; i < NUM_PARS; i++) {
    local_values.vector[i] = values[2+i];
//if (q_index==0) printf("p%d = %g\n",i, local_values.vector[i]);
  }
//if (q_index==0) printf("NUM_VALUES:%d  NUM_PARS:%d  MAX_PD:%d\n", NUM_VALUES, NUM_PARS, MAX_PD);
//if (q_index==0) printf("start:%d stop:%d\n", pd_start, pd_stop);

  double pd_norm = (pd_start == 0 ? 0.0 : result[nq]);
  double this_result = (pd_start == 0 ? 0.0 : result[q_index]);
//if (q_index==0) printf("start %d %g %g\n", pd_start, pd_norm, this_result);

#if MAX_PD>0
  global const double *pd_value = values + NUM_VALUES;
  global const double *pd_weight = pd_value + details->num_weights;
#endif

  // Jump into the middle of the polydispersity loop
#if MAX_PD>4
  int n4=details->pd_length[4];
  int i4=(pd_start/details->pd_stride[4])%n4;
  const int p4=details->pd_par[4];
  global const double *v4 = pd_value + details->pd_offset[4];
  global const double *w4 = pd_weight + details->pd_offset[4];
#endif
#if MAX_PD>3
  int n3=details->pd_length[3];
  int i3=(pd_start/details->pd_stride[3])%n3;
  const int p3=details->pd_par[3];
  global const double *v3 = pd_value + details->pd_offset[3];
  global const double *w3 = pd_weight + details->pd_offset[3];
//if (q_index==0) printf("offset %d: %d %d\n", 3, details->pd_offset[3], NUM_VALUES);
#endif
#if MAX_PD>2
  int n2=details->pd_length[2];
  int i2=(pd_start/details->pd_stride[2])%n2;
  const int p2=details->pd_par[2];
  global const double *v2 = pd_value + details->pd_offset[2];
  global const double *w2 = pd_weight + details->pd_offset[2];
#endif
#if MAX_PD>1
  int n1=details->pd_length[1];
  int i1=(pd_start/details->pd_stride[1])%n1;
  const int p1=details->pd_par[1];
  global const double *v1 = pd_value + details->pd_offset[1];
  global const double *w1 = pd_weight + details->pd_offset[1];
#endif
#if MAX_PD>0
  int n0=details->pd_length[0];
  int i0=(pd_start/details->pd_stride[0])%n0;
  const int p0=details->pd_par[0];
  global const double *v0 = pd_value + details->pd_offset[0];
  global const double *w0 = pd_weight + details->pd_offset[0];
#endif


#if MAX_PD>0
  const int theta_par = details->theta_par;
  const bool fast_theta = (theta_par == p0);
  const bool slow_theta = (theta_par >= 0 && !fast_theta);
  double spherical_correction = 1.0;
#else
  // Note: if not polydisperse the weights cancel and we don't need the
  // spherical correction.
  const double spherical_correction = 1.0;
#endif

  int step = pd_start;


#if MAX_PD>4
  const double weight5 = 1.0;
  while (i4 < n4) {
    local_values.vector[p4] = v4[i4];
    double weight4 = w4[i4] * weight5;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 4, p4, i4, n4, local_values.vector[p4], weight4);
#elif MAX_PD>3
    const double weight4 = 1.0;
#endif
#if MAX_PD>3
  while (i3 < n3) {
    local_values.vector[p3] = v3[i3];
    double weight3 = w3[i3] * weight4;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 3, p3, i3, n3, local_values.vector[p3], weight3);
#elif MAX_PD>2
    const double weight3 = 1.0;
#endif
#if MAX_PD>2
  while (i2 < n2) {
    local_values.vector[p2] = v2[i2];
    double weight2 = w2[i2] * weight3;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 2, p2, i2, n2, local_values.vector[p2], weight2);
#elif MAX_PD>1
    const double weight2 = 1.0;
#endif
#if MAX_PD>1
  while (i1 < n1) {
    local_values.vector[p1] = v1[i1];
    double weight1 = w1[i1] * weight2;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 1, p1, i1, n1, local_values.vector[p1], weight1);
#elif MAX_PD>0
    const double weight1 = 1.0;
#endif
#if MAX_PD>0
  if (slow_theta) { // Theta is not in inner loop
    spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[theta_par])), 1.e-6);
  }
  while(i0 < n0) {
    local_values.vector[p0] = v0[i0];
    double weight0 = w0[i0] * weight1;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 0, p0, i0, n0, local_values.vector[p0], weight0);
    if (fast_theta) { // Theta is in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[p0])), 1.e-6);
    }
#else
    const double weight0 = 1.0;
#endif

//if (q_index == 0) {printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NUM_PARS; i++) printf("p%d=%g ",i, local_values.vector[i]); printf("\n"); }
//if (q_index == 0) printf("sphcor: %g\n", spherical_correction);

    #ifdef INVALID
    if (!INVALID(local_values.table))
    #endif
    {
      // Accumulate I(q)
      // Note: weight==0 must always be excluded
      if (weight0 > cutoff) {
        // spherical correction is set at a minimum of 1e-6, otherwise there
        // would be problems looking at models with theta=90.
        const double weight = weight0 * spherical_correction;
        pd_norm += weight * CALL_VOLUME(local_values.table);

#if defined(MAGNETIC) && NUM_MAGNETIC > 0
        const double qx = q[2*q_index];
        const double qy = q[2*q_index+1];
        const double qsq = qx*qx + qy*qy;

        // Constant across orientation, polydispersity for given qx, qy
        double scattering = 0.0;
        // TODO: what is the magnetic scattering at q=0
        if (qsq > 1.e-16) {
          double p[4];  // dd, du, ud, uu
          p[0] = (qy*cos_mspin + qx*sin_mspin)/qsq;
          p[3] = -p[0];
          p[1] = p[2] = (qy*sin_mspin - qx*cos_mspin)/qsq;

          for (int index=0; index<4; index++) {
            const double xs = spins[index];
            if (xs > 1.e-8) {
              const int spin_flip = (index==1) || (index==2);
              const double pk = p[index];
              for (int axis=0; axis<=spin_flip; axis++) {
                #define M1 NUM_PARS+5
                #define M2 NUM_PARS+8
                #define M3 NUM_PARS+13
                #define SLD(_M_offset, _sld_offset) \
                    local_values.vector[_sld_offset] = xs * (axis \
                    ? (index==1 ? -values[_M_offset+2] : values[_M_offset+2]) \
                    : mag_sld(qx, qy, pk, values[_M_offset], values[_M_offset+1], \
                              (spin_flip ? 0.0 : values[_sld_offset+2])))
                #if NUM_MAGNETIC==1
                    SLD(M1, MAGNETIC_PAR1);
                #elif NUM_MAGNETIC==2
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                #elif NUM_MAGNETIC==3
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                    SLD(M3, MAGNETIC_PAR3);
                #else
                for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                    SLD(M1+3*sk, slds[sk]);
                }
                #endif
                scattering += CALL_IQ(q, q_index, local_values.table);
              }
            }
          }
        }
#else  // !MAGNETIC
//if (q_index==0) printf("calling kernel\n");
        const double scattering = CALL_IQ(q, q_index, local_values.table);
#endif // !MAGNETIC
        this_result += weight * scattering;
      }
    }
    ++step;
#if MAX_PD>0
    if (step >= pd_stop) break;
    ++i0;
  }
  i0 = 0;
#endif
#if MAX_PD>1
    if (step >= pd_stop) break;
    ++i1;
  }
  i1 = 0;
#endif
#if MAX_PD>2
    if (step >= pd_stop) break;
    ++i2;
  }
  i2 = 0;
#endif
#if MAX_PD>3
    if (step >= pd_stop) break;
    ++i3;
  }
  i3 = 0;
#endif
#if MAX_PD>4
    if (step >= pd_stop) break;
    ++i4;
  }
  i4 = 0;
#endif

//if (q_index==0) printf("res: %g/%g\n", this_result, pd_norm);
  // Remember the current result and the updated norm.
  result[q_index] = this_result;
  if (q_index == 0) result[nq] = pd_norm;
}

#undef CALL_IQ
#undef KERNEL_NAME
#endif // cylinder_Iqxy

#if 0 // cylinder_Imagnetic
#define KERNEL_NAME cylinder_Imagnetic
#define MAGNETIC 1
#define CALL_IQ(_q,_i,_v) Iqxy(_q[2*_i],_q[2*_i+1],_v.sld,_v.sld_solvent,_v.radius,_v.length,_v.theta,_v.phi)
#line 1 "/home/pkienzle/src/sasmodels/sasmodels/kernel_iq.cl Imagnetic"

/*
    ##########################################################
    #                                                        #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                              !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNELPY.PY  !!   #
    #   !!                                              !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                        #
    ##########################################################
*/

#ifndef _PAR_BLOCK_ // protected block so we can include this code twice.
#define _PAR_BLOCK_

typedef struct {
#if MAX_PD > 0
    int32_t pd_par[MAX_PD];     // id of the nth polydispersity variable
    int32_t pd_length[MAX_PD];  // length of the nth polydispersity weight vector
    int32_t pd_offset[MAX_PD];  // offset of pd weights in the value & weight vector
    int32_t pd_stride[MAX_PD];  // stride to move to the next index at this level
#endif // MAX_PD > 0
    int32_t num_eval;           // total number of voxels in hypercube
    int32_t num_weights;        // total length of the weights vector
    int32_t num_active;         // number of non-trivial pd loops
    int32_t theta_par;          // id of spherical correction variable
} ProblemDetails;

// Intel HD 4000 needs private arrays to be a multiple of 4 long
typedef struct {
    PARAMETER_TABLE
} ParameterTable;
typedef union {
    ParameterTable table;
    double vector[4*((NUM_PARS+3)/4)];
} ParameterBlock;
#endif // _PAR_BLOCK_


#if defined(MAGNETIC) && NUM_MAGNETIC>0

// Return value restricted between low and high
static double clip(double value, double low, double high)
{
  return (value < low ? low : (value > high ? high : value));
}

// Compute spin cross sections given in_spin and out_spin
// To convert spin cross sections to sld b:
//     uu * (sld - m_sigma_x);
//     dd * (sld + m_sigma_x);
//     ud * (m_sigma_y + 1j*m_sigma_z);
//     du * (m_sigma_y - 1j*m_sigma_z);
static void set_spins(double in_spin, double out_spin, double spins[4])
{
  in_spin = clip(in_spin, 0.0, 1.0);
  out_spin = clip(out_spin, 0.0, 1.0);
  spins[0] = sqrt(sqrt((1.0-in_spin) * (1.0-out_spin))); // dd
  spins[1] = sqrt(sqrt((1.0-in_spin) * out_spin));       // du
  spins[2] = sqrt(sqrt(in_spin * (1.0-out_spin)));       // ud
  spins[3] = sqrt(sqrt(in_spin * out_spin));             // uu
}

static double mag_sld(double qx, double qy, double p,
                       double mx, double my, double sld)
{
    const double perp = qy*mx - qx*my;
    return sld + perp*p;
}

#endif // MAGNETIC

kernel
void KERNEL_NAME(
    int32_t nq,                 // number of q values
    const int32_t pd_start,     // where we are in the polydispersity loop
    const int32_t pd_stop,      // where we are stopping in the polydispersity loop
    global const ProblemDetails *details,
    global const double *values,
    global const double *q, // nq q values, with padding to boundary
    global double *result,  // nq+1 return values, again with padding
    const double cutoff     // cutoff in the polydispersity weight product
    )
{

  // who we are and what element we are working with
  const int q_index = get_global_id(0);
  if (q_index >= nq) return;

  // Storage for the current parameter values.  These will be updated as we
  // walk the polydispersity cube.
  ParameterBlock local_values;

#if defined(MAGNETIC) && NUM_MAGNETIC>0
  // Location of the sld parameters in the parameter vector.
  // These parameters are updated with the effective sld due to magnetism.
  #if NUM_MAGNETIC > 3
  const int32_t slds[] = { MAGNETIC_PARS };
  #endif

  // TODO: could precompute these outside of the kernel.
  // Interpret polarization cross section.
  //     up_frac_i = values[NUM_PARS+2];
  //     up_frac_f = values[NUM_PARS+3];
  //     up_angle = values[NUM_PARS+4];
  double spins[4];
  double cos_mspin, sin_mspin;
  set_spins(values[NUM_PARS+2], values[NUM_PARS+3], spins);
  SINCOS(-values[NUM_PARS+4]*M_PI_180, sin_mspin, cos_mspin);
#endif // MAGNETIC

  // Fill in the initial variables
  //   values[0] is scale
  //   values[1] is background
  for (int i=0; i < NUM_PARS; i++) {
    local_values.vector[i] = values[2+i];
//if (q_index==0) printf("p%d = %g\n",i, local_values.vector[i]);
  }
//if (q_index==0) printf("NUM_VALUES:%d  NUM_PARS:%d  MAX_PD:%d\n", NUM_VALUES, NUM_PARS, MAX_PD);
//if (q_index==0) printf("start:%d stop:%d\n", pd_start, pd_stop);

  double pd_norm = (pd_start == 0 ? 0.0 : result[nq]);
  double this_result = (pd_start == 0 ? 0.0 : result[q_index]);
//if (q_index==0) printf("start %d %g %g\n", pd_start, pd_norm, this_result);

#if MAX_PD>0
  global const double *pd_value = values + NUM_VALUES;
  global const double *pd_weight = pd_value + details->num_weights;
#endif

  // Jump into the middle of the polydispersity loop
#if MAX_PD>4
  int n4=details->pd_length[4];
  int i4=(pd_start/details->pd_stride[4])%n4;
  const int p4=details->pd_par[4];
  global const double *v4 = pd_value + details->pd_offset[4];
  global const double *w4 = pd_weight + details->pd_offset[4];
#endif
#if MAX_PD>3
  int n3=details->pd_length[3];
  int i3=(pd_start/details->pd_stride[3])%n3;
  const int p3=details->pd_par[3];
  global const double *v3 = pd_value + details->pd_offset[3];
  global const double *w3 = pd_weight + details->pd_offset[3];
//if (q_index==0) printf("offset %d: %d %d\n", 3, details->pd_offset[3], NUM_VALUES);
#endif
#if MAX_PD>2
  int n2=details->pd_length[2];
  int i2=(pd_start/details->pd_stride[2])%n2;
  const int p2=details->pd_par[2];
  global const double *v2 = pd_value + details->pd_offset[2];
  global const double *w2 = pd_weight + details->pd_offset[2];
#endif
#if MAX_PD>1
  int n1=details->pd_length[1];
  int i1=(pd_start/details->pd_stride[1])%n1;
  const int p1=details->pd_par[1];
  global const double *v1 = pd_value + details->pd_offset[1];
  global const double *w1 = pd_weight + details->pd_offset[1];
#endif
#if MAX_PD>0
  int n0=details->pd_length[0];
  int i0=(pd_start/details->pd_stride[0])%n0;
  const int p0=details->pd_par[0];
  global const double *v0 = pd_value + details->pd_offset[0];
  global const double *w0 = pd_weight + details->pd_offset[0];
#endif


#if MAX_PD>0
  const int theta_par = details->theta_par;
  const bool fast_theta = (theta_par == p0);
  const bool slow_theta = (theta_par >= 0 && !fast_theta);
  double spherical_correction = 1.0;
#else
  // Note: if not polydisperse the weights cancel and we don't need the
  // spherical correction.
  const double spherical_correction = 1.0;
#endif

  int step = pd_start;


#if MAX_PD>4
  const double weight5 = 1.0;
  while (i4 < n4) {
    local_values.vector[p4] = v4[i4];
    double weight4 = w4[i4] * weight5;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 4, p4, i4, n4, local_values.vector[p4], weight4);
#elif MAX_PD>3
    const double weight4 = 1.0;
#endif
#if MAX_PD>3
  while (i3 < n3) {
    local_values.vector[p3] = v3[i3];
    double weight3 = w3[i3] * weight4;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 3, p3, i3, n3, local_values.vector[p3], weight3);
#elif MAX_PD>2
    const double weight3 = 1.0;
#endif
#if MAX_PD>2
  while (i2 < n2) {
    local_values.vector[p2] = v2[i2];
    double weight2 = w2[i2] * weight3;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 2, p2, i2, n2, local_values.vector[p2], weight2);
#elif MAX_PD>1
    const double weight2 = 1.0;
#endif
#if MAX_PD>1
  while (i1 < n1) {
    local_values.vector[p1] = v1[i1];
    double weight1 = w1[i1] * weight2;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 1, p1, i1, n1, local_values.vector[p1], weight1);
#elif MAX_PD>0
    const double weight1 = 1.0;
#endif
#if MAX_PD>0
  if (slow_theta) { // Theta is not in inner loop
    spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[theta_par])), 1.e-6);
  }
  while(i0 < n0) {
    local_values.vector[p0] = v0[i0];
    double weight0 = w0[i0] * weight1;
//if (q_index == 0) printf("step:%d level %d: p:%d i:%d n:%d value:%g weight:%g\n", step, 0, p0, i0, n0, local_values.vector[p0], weight0);
    if (fast_theta) { // Theta is in inner loop
      spherical_correction = fmax(fabs(cos(M_PI_180*local_values.vector[p0])), 1.e-6);
    }
#else
    const double weight0 = 1.0;
#endif

//if (q_index == 0) {printf("step:%d of %d, pars:",step,pd_stop); for (int i=0; i < NUM_PARS; i++) printf("p%d=%g ",i, local_values.vector[i]); printf("\n"); }
//if (q_index == 0) printf("sphcor: %g\n", spherical_correction);

    #ifdef INVALID
    if (!INVALID(local_values.table))
    #endif
    {
      // Accumulate I(q)
      // Note: weight==0 must always be excluded
      if (weight0 > cutoff) {
        // spherical correction is set at a minimum of 1e-6, otherwise there
        // would be problems looking at models with theta=90.
        const double weight = weight0 * spherical_correction;
        pd_norm += weight * CALL_VOLUME(local_values.table);

#if defined(MAGNETIC) && NUM_MAGNETIC > 0
        const double qx = q[2*q_index];
        const double qy = q[2*q_index+1];
        const double qsq = qx*qx + qy*qy;

        // Constant across orientation, polydispersity for given qx, qy
        double scattering = 0.0;
        // TODO: what is the magnetic scattering at q=0
        if (qsq > 1.e-16) {
          double p[4];  // dd, du, ud, uu
          p[0] = (qy*cos_mspin + qx*sin_mspin)/qsq;
          p[3] = -p[0];
          p[1] = p[2] = (qy*sin_mspin - qx*cos_mspin)/qsq;

          for (int index=0; index<4; index++) {
            const double xs = spins[index];
            if (xs > 1.e-8) {
              const int spin_flip = (index==1) || (index==2);
              const double pk = p[index];
              for (int axis=0; axis<=spin_flip; axis++) {
                #define M1 NUM_PARS+5
                #define M2 NUM_PARS+8
                #define M3 NUM_PARS+13
                #define SLD(_M_offset, _sld_offset) \
                    local_values.vector[_sld_offset] = xs * (axis \
                    ? (index==1 ? -values[_M_offset+2] : values[_M_offset+2]) \
                    : mag_sld(qx, qy, pk, values[_M_offset], values[_M_offset+1], \
                              (spin_flip ? 0.0 : values[_sld_offset+2])))
                #if NUM_MAGNETIC==1
                    SLD(M1, MAGNETIC_PAR1);
                #elif NUM_MAGNETIC==2
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                #elif NUM_MAGNETIC==3
                    SLD(M1, MAGNETIC_PAR1);
                    SLD(M2, MAGNETIC_PAR2);
                    SLD(M3, MAGNETIC_PAR3);
                #else
                for (int sk=0; sk<NUM_MAGNETIC; sk++) {
                    SLD(M1+3*sk, slds[sk]);
                }
                #endif
                scattering += CALL_IQ(q, q_index, local_values.table);
              }
            }
          }
        }
#else  // !MAGNETIC
        const double scattering = CALL_IQ(q, q_index, local_values.table);
#endif // !MAGNETIC
        this_result += weight * scattering;
      }
    }
    ++step;
#if MAX_PD>0
    if (step >= pd_stop) break;
    ++i0;
  }
  i0 = 0;
#endif
#if MAX_PD>1
    if (step >= pd_stop) break;
    ++i1;
  }
  i1 = 0;
#endif
#if MAX_PD>2
    if (step >= pd_stop) break;
    ++i2;
  }
  i2 = 0;
#endif
#if MAX_PD>3
    if (step >= pd_stop) break;
    ++i3;
  }
  i3 = 0;
#endif
#if MAX_PD>4
    if (step >= pd_stop) break;
    ++i4;
  }
  i4 = 0;
#endif

//if (q_index==0) printf("res: %g/%g\n", this_result, pd_norm);
  // Remember the current result and the updated norm.
  result[q_index] = this_result;
  if (q_index == 0) result[nq] = pd_norm;
}

#undef MAGNETIC
#undef CALL_IQ
#undef KERNEL_NAME
#endif // cylinder_magnetic
