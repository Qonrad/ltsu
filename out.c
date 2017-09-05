/* out.c */
#include <math.h>
#include <stdio.h>
#include "lts.h"

#define T_RES 1.0
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
static double t_old,yold,yolder,yprev,tprev;
static double yprev2,yold2,yolder2;
extern double current[C];
static double oldcurrent[C];
static double oldstate[N];
static double last_t=0,last_t2;
double tcross,tcross1,tcross2,period[2];
FILE *pi;

int i;
if(*nr==1) 
{
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
    pi= fopen("interval.data","w");
    fclose(pi);
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
if(((yold < THRESHOLD) && (state[0] > THRESHOLD) ) && ((yold2 > THRESHOLD) || (state[6] < THRESHOLD))) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
    period[0] = tcross - last_t;
    pi= fopen("interval.data","a");
    fprintf(pi,"%d %f\n",0,period[0]);
    fclose(pi);
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
    yprev=state[0];
    tprev=*x;
    last_t=tcross;}
if( ((yold2 < THRESHOLD) && (state[6] > THRESHOLD)) && ((yold > THRESHOLD) || (state[0] < THRESHOLD)) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold2)/(state[6] - yold2);
    period[1] = tcross - last_t;
    pi= fopen("interval.data","a");
    fprintf(pi,"%d %f\n",1,period[1]);
    fclose(pi);
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
    yprev2=state[6];
    tprev=*x;
    last_t=tcross;}
if( ((yold2 < THRESHOLD) && (state[6] > THRESHOLD)) && ((yold < THRESHOLD) && (state[0] > THRESHOLD)) ) 
   { tcross2 = t_old + (*x - t_old)*(THRESHOLD - yold2)/(state[6] - yold2);
     tcross1 = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
      if(tcross1 <= tcross2)
      {
    period[0] = tcross1 - last_t;
    period[1] = tcross2 - tcross1;
    pi= fopen("interval.data","a");
    fprintf(pi,"%d %f\n",0,period[0]);
    fprintf(pi,"%d %f\n",1,period[1]);
    fclose(pi);
    last_t = tcross2;
      }
      if(tcross1 >  tcross2)
      {
    period[1] = tcross2 - last_t;
    period[0] = tcross1 - tcross2;
    pi= fopen("interval.data","a");
    fprintf(pi,"%d %f\n",1,period[1]);
    fprintf(pi,"%d %f\n",0,period[0]);
    fclose(pi);
    last_t = tcross1;
      }
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
    yprev=state[0];
    yprev2=state[6];
    tprev=*x;
    }
if(*nr!=1 && t_old> tprev && (((t_old - tprev) > T_RES) || (yold-yolder)/(state[0]-yold)<0
|| (yold2-yolder2)/(state[6]-yold2)<0
|| fabs(yold-yprev)> V_LOW_RES || fabs(yold2-yprev2)> V_LOW_RES ))
{
     printf("%f ",t_old);
     for(i=0;i<N;i++) printf("%f ", oldstate[i]);
     printf("\n");
 fflush(stdout); 
yprev=yold;
yprev2=yold2;
tprev=t_old;}
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<N;i++) oldcurrent[i]=current[i];
}
