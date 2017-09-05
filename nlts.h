#define THRESHOLD -14.0 /* -14.0 */
#define DELAY 30.48
#define INCREMENT 100 
#define N 6
#define NN 20
#define C 5
#define V1 0  /* mV */
#define M1 1  /* mV */
#define H1 2  /* mM */
#define N1 3 
#define MN1 4 
#define S1 5 
#define V2 6 
#define M2 7 
#define H2 8 
#define N2 9 
#define MN2 10 
#define S2 11 
#define I_NA1 0  
#define I_K1 1 
#define I_M1 2 
#define I_L1 3
#define I_S1 4 
#define I_NA2 5
#define I_K2 6 
#define I_M2 7 
#define I_L2 8 
#define I_S2 9 
#define CM 1.00 /*uF/cm2*/
#define I_APP 2.0 /*uA/cm2*/
#define EPSILON 0.00  /* 0.07 0.14 0.2 */
#define E_NA  50.0
#define E_K  -100.0
#define E_L  -67.0
#define E_SYN  -80.0
#define G_NA 100.0 /* mS/cm2*/
#define G_K   80.0
#define G_M   2.0
#define G_L   0.1
#define G_SYN  0.25/NN /*0.25/9.0  0.35 0.250 mS/cm2*/
#define TAUSYN 10.0 /*  1.0 3.0 2.0   */


/* INTEGRATION PARAMETERS */

/* these only apply to C code */
#define LWORK 2*58567 /* 679  N(4*N + 8)+7*/
#define LIWORK 2*367  /*43 3*N+7*/
#define LRCONT  2*484 /* 52 4*N+4*/
#define DURATION 5
#define TIMEON 0
#define TIMEOFF 0
#define PERIOD  250
#define ENDTIME 5000
#define START_TIME  0

#define PRINT  0
#define DEBUG 1
#define PULSEON 300000
#define PULSEOFF 300000


