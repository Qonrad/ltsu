/* fix epsilon!*/
#include <math.h>
#include <stdio.h>
#include "nlts.h"




double boltz(v,half,slope)
double v,half,slope;
{
double arg;
arg = -(v-half)/slope;
if ((arg>50.0) || (arg<-50.0))
{if (arg>50.0) return(0.0);
else  return(1.0);}
else return(1.0/(1.0 + exp(arg)));}



int deriv_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
extern double current[C*NN]; 
int el,n,j;
double time_,epsilon;
time_ = *xp;

el = *np;

for(j=0;j<NN;j++)
{
current[I_NA1 + C*j] = G_NA*Y[H1 + N*j]*pow(Y[M1 + N*j],3.0)*(Y[V1 + N*j]-E_NA);
current[I_K1 + C*j] = G_K*pow(Y[N1 + N*j],4.0)*(Y[V1 + N*j] - E_K);
current[I_M1 + C*j] = G_M*Y[MN1 + N*j]*(Y[V1 + N*j] - E_K);
current[I_L1 + C*j] =  G_L*(Y[V1 + N*j] - E_L);
current[I_S1 + C*j] = 0.0;
for (n=0;n<NN;n++)
 {
//if(time_<60&n!=j) current[I_S1 + C*j] =   current[I_S1 + C*j] + G_SYN*Y[S1 +N*n]*(Y[V1 +N*j] - E_SYN);
if(n!=j) current[I_S1 + C*j] =   current[I_S1 + C*j] + G_SYN*Y[S1 +N*n]*(Y[V1 +N*j] - E_SYN);
 }
F[V1 + N*j] = (I_APP + EPSILON - current[I_NA1 + C*j] - current[I_K1 + C*j] - current[I_M1 + C*j] - current[I_L1 + C*j] - current[I_S1 + C*j])/CM;
F[M1 + N*j] = 0.32*(Y[V1 + N*j]+54.0)/(1.0-exp(-(Y[V1 + N*j]+54.0)/4.0))*(1.0-Y[M1 + N*j])-0.28*(Y[V1 + N*j]+27.0)/(exp((Y[V1 + N*j]+27.0)/5.0)-1.0)*Y[M1 + N*j];   
F[H1 + N*j] = 0.128*exp(-(Y[V1 + N*j]+50.0)/18.0)*(1.0-Y[H1 + N*j])-4.0/(1.0+exp(-(Y[V1 + N*j]+27.0)/5.0))*Y[H1 + N*j];   
F[N1 + N*j] = 0.032*(Y[V1 + N*j]+52.0)/(1.0-exp(-(Y[V1 + N*j]+52.0)/5.0))*(1.0-Y[N1 + N*j])-0.5*exp(-(Y[V1 + N*j]+57.0)/40.0)*Y[N1 + N*j];  
F[MN1 + N*j] = 3.209*0.0001*((Y[V1 + N*j]+30.0)/(1.0-exp(-(Y[V1 + N*j]+30.0)/9.0))*(1.0-Y[MN1 + N*j]) 
              + (Y[V1 + N*j]+30.0)/(1.0-exp((Y[V1 + N*j]+30.0)/9.0))*Y[MN1 + N*j]); 
F[S1 + N*j] = 2*(1+tanh(Y[V1 + N*j]/4.0))*(1-Y[S1 + N*j])-Y[S1 + N*j]/TAUSYN;  
 }

return 0;
   }


void scan_(Y) 
double Y[N*NN];
{FILE *fopen(),*sp;
int i;
sp = fopen("state.data","r");
for(i=0;i<N*NN;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);}

void dump_(Y) 
double Y[N*NN];
{FILE *fopen(),*sp;
int i;
sp = fopen("end.data","w");
for(i=0;i<N*NN;i++) fprintf(sp,"%.16f\n",Y[i]);
fclose(sp);}

int mas(n,amas,l)
        int *n;
        double *amas;
        int *l;
{return 0;}

int dummy(n,t,y,ydot)
        int *n;
        double *t;
        double *y;
        double *ydot;
{return 0;}

