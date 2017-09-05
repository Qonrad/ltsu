#include <stdio.h>
#include <math.h>
#include "nlts.h"
#include "nlts.c"
#include "cblock.h"
#include <malloc.h>


/* global variable declarations */

 
double current[C*NN],vset1;
 
void main()
{
double state[N*NN];
double fstate[N*NN];
double rtol[N*NN],atol[N*NN];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int *iwork;
double *work;
extern int out_();
int n,idid;

double x,xend,h;
FILE *fp;
int i,j;
n=N*NN;
for(i=0;i<N*NN;i++) {rtol[i]=1.0e-6; //e-8
                  atol[i]=1.0e-9;} //e-11

for(j=0;j<NN;j++) {rtol[j*N]=1.0e-8; //e-8
                  atol[j*N]=1.0e-11;} //e-11
itol=1;
ijac=0;
mljac=n;
mujac=0;
imas=0;
mlmas=n;
mumas=0;
lwork= n*(4*n +8) +7;
liwork=3*n +7;
lrcont=4*n + 4;
iout=1;

work = malloc((lwork) * sizeof(double));
iwork = malloc((liwork) * sizeof(int));
for(i=0;i<7;i++) {
	iwork[i]=0;
    work[i]=0.0;
   }

iwork[2] = 20000000; //20000000
h=1e-4; // 1e-4
x=START_TIME;
scan_(state);
deriv_(&n,&x,state,fstate);
if (PULSEON>0 && PULSEON<ENDTIME) {xend = PULSEON;
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
//Where is the loop holding this all together?, Why doesn't it all just run once?
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

