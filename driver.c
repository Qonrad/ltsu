#include <stdio.h>
#include <math.h>
#include "lts.h"
#include "lts.c"
#include "cblock.h"


/* global variable declarations */

/* double F[N];
double Y[N];
long int *np;
double *xp; */
 
double current[C],vset1;
 
void main()
{
double state[N];
double fstate[N];
double rtol[N],atol[N];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int iwork[LIWORK];
double work[LWORK];
extern int out_();
int n,idid;

double x,xend,h;
FILE *fp;
int i,j;
n=N;
for(i=0;i<7;i++) {iwork[i]=0;
                  work[i]=0.0;}
for(i=1;i<N;i++) {rtol[i]=1.0e-6;
                  atol[i]=1.0e-11;}
rtol[0]=1.0e-8; 
atol[0]=1.0e-14;
rtol[6]=1.0e-8; 
atol[6]=1.0e-14;  

iwork[2] = 10000000;
itol=1;
ijac=0;
mljac=n;
mujac=0;
imas=0;
mlmas=n;
mumas=0;
lwork=LWORK;
liwork=LIWORK;
lrcont=LRCONT;
iout=1;

h=1e-4;
x=START_TIME;
scan_(state);
deriv_(&n,&x,state,fstate);
if (PULSEON>0 && PULSEON<ENDTIME)
{xend = PULSEON;
radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
}

if (xend<PULSEOFF&& PULSEOFF<ENDTIME)
{xend = PULSEOFF;
radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
}

if(x<ENDTIME)
{xend=ENDTIME;

radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);}

dump_(state);
fflush(stdout);
}

