/* out.c */
#include <math.h>
#include <stdio.h>
#include "lts.h"

#define T_RES 1
#define V_LOW_RES 2.0
#define V_HIGH_RES 0.05
#define SAMPLE 4

int out_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
static double oldcurrent[C];
static double oldstate[N];
extern double phi();
int i;
if(*nr==1) 
{
/*    printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");*/
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
/*if( (yold > yolder) && (state[0] < yold) && (yold > 0)) 
   {printf("%f %f\n",t_old,t_old - last_t);
    last_t = t_old;}*/
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
/* printf("%f %f\n",tcross,tcross - last_t);*/
   period[0] = tcross - last_t;
   initial[0] = THRESHOLD;
   initial[1] = oldstate[1] + (state[1] - oldstate[1])*(THRESHOLD - yold)/(state[0] - yold);
   initial[2] = oldstate[2] + (state[2] - oldstate[2])*(THRESHOLD - yold)/(state[0] - yold);
   initial[3] = oldstate[3] + (state[3] - oldstate[3])*(THRESHOLD - yold)/(state[0] - yold);
   initial[4] = oldstate[4] + (state[4] - oldstate[4])*(THRESHOLD - yold)/(state[0] - yold);
   initial[5] = 0.000000;
    last_t = tcross;}
if( (yold2 < THRESHOLD) && (state[6] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold2)/(state[6] - yold2);
   /*printf("%f %f\n",tcross,tcross - last_t2);*/
   period[1] = tcross - last_t2;
   initial[6] = THRESHOLD;
   initial[7] = oldstate[7] + (state[7] - oldstate[7])*(THRESHOLD - yold2)/(state[6] - yold2);
   initial[8] = oldstate[8] + (state[8] - oldstate[8])*(THRESHOLD - yold2)/(state[6] - yold2);
   initial[9] = oldstate[9] + (state[9] - oldstate[9])*(THRESHOLD - yold2)/(state[6] - yold2);
   initial[10] = oldstate[10] + (state[10] - oldstate[10])*(THRESHOLD - yold2)/(state[6] - yold2);
   initial[11] =  0.000000;
    last_t2 = tcross;}
/*if( (yold2 > yolder2)&&(state[6] < yold2) && (yold2 > 0))
   {printf("%f %f\n",t_old,t_old - last_t2);
    last_t2 =  t_old;}*/
if((fabs(tprev - t_old) > T_RES) || (yold-yolder)/(state[0]-yold)<0
|| (yold2-yolder2)/(state[6]-yold2)<0
|| fabs(yold-yprev)> V_LOW_RES || fabs(yold2-yprev2)> V_LOW_RES )
{
    /* printf("%f ",t_old);
     for(i=0;i<N;i++) printf("%f ", oldstate[i]);
     printf("\n");
 fflush(stdout); */
yprev=yold;
yprev2=yold2;
tprev=t_old;}
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}
int altout_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
extern double f1,f2,f3;
extern int flag1,flag2,flag3,flag4;
static double oldcurrent[C];
static double oldstate[N];
extern double phi();
int i;
if(*nr==1) 
{
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   if((flag1==1)&& (flag2==0)) f1 = (tcross - period[0])/period[0];
   if((flag2==1)&& (flag3==0)) f2 = (tcross - last_t - period[0])/period[0];
   if((flag3==1)&& (!flag4)) f3 = (tcross -last_t - period[0])/period[0];
    last_t=tcross;}
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}

int altout2_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
extern double f1,f2,f3;
extern int flag1,flag2,flag3,flag4;
static double oldcurrent[C];
static double oldstate[N];
extern double phi();
int i;
if(*nr==1) 
{
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
if( (yold2 < THRESHOLD) && (state[6] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold2)/(state[6] - yold2);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   if((flag1==1)&& (flag2==0)) f1 = (tcross - period[1])/period[1];
   if((flag2==1)&& (flag3==0)) f2 = (tcross - last_t - period[1])/period[1];
   if((flag3==1)&&(!flag4)) f3 = (tcross -last_t - period[1])/period[1];
    last_t=tcross;}
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}


int oldout_(nr,xold,x,state,n,irtrn)
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
extern double phi();
int i;
if(*nr==1) 
{
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
else if((fabs(tprev - t_old) > T_RES) || (yold-yolder)/(state[0]-yold)<0
|| (yold2-yolder2)/(state[6]-yold2)<0
|| fabs(yold-yprev)> V_LOW_RES || fabs(yold2-yprev2)> V_LOW_RES )
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
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}

int debugout_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
static double oldcurrent[C];
static double oldstate[N];
extern int flag1, flag2, flag3,flag4;
extern double f1, f2, f3;
extern double phi();
int i;
if(*nr==1) 
{
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
else if((fabs(tprev - t_old) > T_RES) || (yold-yolder)/(state[0]-yold)<0
|| (yold2-yolder2)/(state[6]-yold2)<0
|| fabs(yold-yprev)> V_LOW_RES || fabs(yold2-yprev2)> V_LOW_RES )
{
     printf("%f ",t_old);
     for(i=0;i<N;i++) printf("%f ", oldstate[i]);
     printf("\n");
 fflush(stdout); 
yprev=yold;
yprev2=yold2;
tprev=t_old;}
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   if ((flag1==1)&& (flag2==0)) f1 = (tcross - period[0])/period[0];
   if ((flag2==1)&& (flag3==0)) f2 = (tcross - last_t - period[0])/period[0];
   if ((flag3==1) &&(!flag4)) f3 = (tcross -last_t - period[0])/period[0];
    last_t=tcross;}
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}
int debugout2_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],initial[N],period[2];
static double oldcurrent[C];
static double oldstate[N];
extern int flag1, flag2, flag3,flag4;
extern double f1, f2, f3;
extern double phi();
int i;
if(*nr==1) 
{
     printf("%f ",*x);
     for(i=0;i<N;i++) printf("%f ", state[i]);
     printf("\n");
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
yprev=yold;
t_old=*x;}
else if((fabs(tprev - t_old) > T_RES) || (yold-yolder)/(state[0]-yold)<0
|| (yold2-yolder2)/(state[6]-yold2)<0
|| fabs(yold-yprev)> V_LOW_RES || fabs(yold2-yprev2)> V_LOW_RES )
{
     printf("%f ",t_old);
     for(i=0;i<N;i++) printf("%f ", oldstate[i]);
     printf("\n");
 fflush(stdout); 
yprev=yold;
yprev2=yold2;
tprev=t_old;}
if( (yold2 < THRESHOLD) && (state[6] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold2)/(state[6] - yold2);
   if(flag3) flag4=1;
   if(flag2) flag3=1;
   if(flag1) flag2=1;
   flag1 = 1;
   if((flag1==1)&& (flag2==0)) f1 = (tcross - period[1])/period[1];
   if((flag2==1)&& (flag3==0)) f2 = (tcross - last_t - period[1])/period[1];
   if ((flag3==1) && (!flag4)) f3 = (tcross -last_t - period[1])/period[1];
    last_t=tcross;}
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}
int sout_(nr,xold,x,state,n,irtrn)
        int *nr;
        double *xold;
        double *x;
        double *state;
        int *n;
        int *irtrn;
{
static double t_old,yold,yolder,yprev,tprev;
static double tcross,yprev2,yold2,yolder2,last_t,last_t2;
extern double current[C],period[2];
static double oldcurrent[C];
static double oldstate[N];
extern double phi();
static int set = 0;
double sample[N];
double chunk;
int i,j;
if( (yold < THRESHOLD) && (state[0] > THRESHOLD) ) 
   { tcross = t_old + (*x - t_old)*(THRESHOLD - yold)/(state[0] - yold);
   sample[0] = THRESHOLD;
   sample[1] = oldstate[1] + (state[1] - oldstate[1])*(THRESHOLD - yold)/(state[0] - yold);
   sample[2] = oldstate[2] + (state[2] - oldstate[2])*(THRESHOLD - yold)/(state[0] - yold);
   sample[3] = oldstate[3] + (state[3] - oldstate[3])*(THRESHOLD - yold)/(state[0] - yold);
   sample[4] = oldstate[4] + (state[4] - oldstate[4])*(THRESHOLD - yold)/(state[0] - yold);
   sample[5] = 0.000000;
   for(j=0;j<4;j++) printf("%f ",sample[j]);
   printf("\n");
   set=1;
   }
chunk = 1.0/(1.0*SAMPLE);
for(i=1;i<SAMPLE+1;i++)if( set && (*x-tcross >= chunk*i*period[0] && t_old - tcross < chunk*i*period[0]))
   { 
   sample[0] = oldstate[0] + (state[0] - oldstate[0])*(chunk*i*period[0] - t_old + tcross)/(*x - t_old);
   sample[1] = oldstate[1] + (state[1] - oldstate[1])*(chunk*i*period[0] - t_old + tcross)/(*x - t_old);
   sample[2] = oldstate[2] + (state[2] - oldstate[2])*(chunk*i*period[0] - t_old + tcross)/(*x - t_old);
   sample[3] = oldstate[3] + (state[3] - oldstate[3])*(chunk*i*period[0] - t_old + tcross)/(*x - t_old);
   sample[4] = oldstate[4] + (state[4] - oldstate[4])*(chunk*i*period[0] - t_old + tcross)/(*x - t_old);
   sample[5] = 0.000000;
   for(j=0;j<6;j++) printf("%f ",sample[j]);
   printf("\n");
   }
yolder=yold;
yolder2=yold2;
yold=state[0];
yold2=state[6];
t_old=*x;
for(i=0;i<N;i++) oldstate[i]=state[i];
for(i=0;i<C;i++) oldcurrent[i]=current[i];
}
