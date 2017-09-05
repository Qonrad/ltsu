/* out.c */
#include <math.h>
#include <stdio.h>
#include "nlts.h"

#define T_RES 0.1
#define V_LOW_RES 2.0
#define V_HIGH_RES 0.05

int out_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold[NN],yolder[NN],yprev[NN],tprev,lastt=0.0;
extern double current[C*NN];
static double oldcurrent[C*NN];
static double oldstate[N*NN];
static double last_t[NN];
double tcross[NN],period[NN],interval;
FILE *pi;

int i,j;
if(*nr==1) 
{
     printf("%f ",*x);
     for(i=0;i<N*NN;i++) printf("%f ", state[i]);
     //for(i=0;i<C*NN;i++) printf("%f ", current[i]);
     //for(i=0;i<NN;i++) printf("%f ", current[4+i*C]);
     printf("\n");
    pi= fopen("interval.data","w");
    fclose(pi);
for(i=0;i<NN;i++)
{
yolder[i]=yold[i];
yold[i]=state[N*i];
yprev[i]=yold[i];
}
t_old=*x;}
/* what if there are multiple crossings in a loop? */
 for(j=0;j<NN;j++)
 {
if((yold[j] < THRESHOLD) && (state[N*j] > THRESHOLD) ) 
   { 
    tcross[j] = t_old + (*x - t_old)*(THRESHOLD - yold[j])/(state[N*j] - yold[j]);
   if(lastt > 0.00001)
   {
    period[j] = tcross[j] - last_t[j];
   interval = tcross[j] - lastt;
    pi= fopen("interval.data","a");
    fprintf(pi,"%d %f\n",j,interval);
    fclose(pi);
     printf("%f ",*x);
     for(i=0;i<N*NN;i++) printf("%f ", state[i]);
     // for(i=0;i<C*NN;i++) printf("%f ", current[i]);
     //for(i=0;i<NN;i++) printf("%f ", current[4+i*C]);
     printf("\n");
    tprev=*x;
   }
last_t[j]= tcross[j];
lastt= tcross[j];
    for(i=0;i<NN;i++) yprev[i]=state[N*i];
    }
    }
if(*nr!=1 && *x> tprev && ((*x - tprev) > T_RES) )
{
     printf("%f ",*x);
     for(i=0;i<N*NN;i++) printf("%f ", state[i]);
     // for(i=0;i<C*NN;i++) printf("%f ", current[i]);
     //for(i=0;i<NN;i++) printf("%f ", current[4+i*C]);
     printf("\n");
 fflush(stdout); 
for(i=0;i<NN;i++)
{
yprev[i]=yold[i];
}
tprev=*x;}
for(i=0;i<NN;i++)
{
yolder[i]=yold[i];
yold[i]=state[N*i];
}
t_old=*x;
for(i=0;i<N*NN;i++) oldstate[i]=state[i];
for(i=0;i<C*NN;i++) oldcurrent[i]=current[i];
}
