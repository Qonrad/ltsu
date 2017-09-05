#include <math.h>
#include <stdio.h>
#include "lts.h"




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
extern double gsyn1,gsyn2,current[C]; 
int el;
double time_;
time_ = *xp;

el = *np;

current[I_NA1] = G_NA*Y[H1]*pow(Y[M1],3.0)*(Y[V1]-E_NA);
current[I_K1] = G_K*pow(Y[N1],4.0)*(Y[V1] - E_K);
current[I_M1] = G_M*Y[MN1]*(Y[V1] - E_K);
current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_S1] =   gsyn1*Y[S2]*(Y[V1] - E_SYN);
F[V1] = (I_APP + EPSILON - current[I_NA1] - current[I_K1] - current[I_M1] - current[I_L1] - current[I_S1])/CM;
F[M1] = 0.32*(Y[V1]+54.0)/(1.0-exp(-(Y[V1]+54.0)/4.0))*(1.0-Y[M1])-0.28*(Y[V1]+27.0)/(exp((Y[V1]+27.0)/5.0)-1.0)*Y[M1];   
F[H1] = 0.128*exp(-(Y[V1]+50.0)/18.0)*(1.0-Y[H1])-4.0/(1.0+exp(-(Y[V1]+27.0)/5.0))*Y[H1];   
F[N1] = 0.032*(Y[V1]+52.0)/(1.0-exp(-(Y[V1]+52.0)/5.0))*(1.0-Y[N1])-0.5*exp(-(Y[V1]+57.0)/40.0)*Y[N1];  
F[MN1] = 3.209*0.0001*((Y[V1]+30.0)/(1.0-exp(-(Y[V1]+30.0)/9.0))*(1.0-Y[MN1]) + (Y[V1]+30.0)/(1.0-exp((Y[V1]+30.0)/9.0))*Y[MN1]); 
F[S1] = 2*(1+tanh(Y[V2]/4.0))*(1-Y[S1])-Y[S1]/TAUSYN;  

current[I_NA2] = G_NA*Y[H2]*pow(Y[M2],3.0)*(Y[V2]-E_NA);
current[I_K2] = G_K*pow(Y[N2],4.0)*(Y[V2] - E_K);
current[I_M2] = G_M*Y[MN2]*(Y[V2] - E_K);
current[I_L2] =  G_L*(Y[V2] - E_L);
current[I_S2] =   gsyn2*Y[S1]*(Y[V2] - E_SYN);
F[V2] = (I_APP + EPSILON - current[I_NA2] - current[I_K2] - current[I_M2] - current[I_L2] - current[I_S2])/CM;
F[M2] = 0.32*(Y[V2]+54.0)/(1.0-exp(-(Y[V2]+54.0)/4.0))*(1.0-Y[M2])-0.28*(Y[V2]+27.0)/(exp((Y[V2]+27.0)/5.0)-1.0)*Y[M2];   
F[H2] = 0.128*exp(-(Y[V2]+50.0)/18.0)*(1.0-Y[H2])-4.0/(1.0+exp(-(Y[V2]+27.0)/5.0))*Y[H2];   
F[N2] = 0.032*(Y[V2]+52.0)/(1.0-exp(-(Y[V2]+52.0)/5.0))*(1.0-Y[N2])-0.5*exp(-(Y[V2]+57.0)/40.0)*Y[N2];  
F[MN2] = 3.209*0.0001*((Y[V2]+30.0)/(1.0-exp(-(Y[V2]+30.0)/9.0))*(1.0-Y[MN2]) + (Y[V2]+30.0)/(1.0-exp((Y[V2]+30.0)/9.0))*Y[MN2]); 
F[S2] = 2*(1+tanh(Y[V2]/4.0))*(1-Y[S2])-Y[S2]/TAUSYN;  
return 0;
   }

int deriv1_(np,xp,Y,F)
double *F,*Y;
int *np;
double *xp;
{
extern double gsyn1,gsyn2,current[C]; 
int el;
double time_;
time_ = *xp;

el = *np;

current[I_NA1] = G_NA*Y[H1]*pow(Y[M1],3.0)*(Y[V1]-E_NA);
current[I_K1] = G_K*pow(Y[N1],4.0)*(Y[V1] - E_K);
current[I_M1] = G_M*Y[MN1]*(Y[V1] - E_K);
current[I_L1] =  G_L*(Y[V1] - E_L);
current[I_S1] =   gsyn1*Y[S2]*(Y[V1] - E_SYN) +G_SYN2*Y[S1]*(Y[V1] - E_SYN);
F[V1] = (I_APP + EPSILON - current[I_NA1] - current[I_K1] - current[I_M1] - current[I_L1] - current[I_S1])/CM;
F[M1] = 0.32*(Y[V1]+54.0)/(1.0-exp(-(Y[V1]+54.0)/4.0))*(1.0-Y[M1])-0.28*(Y[V1]+27.0)/(exp((Y[V1]+27.0)/5.0)-1.0)*Y[M1];   
F[H1] = 0.128*exp(-(Y[V1]+50.0)/18.0)*(1.0-Y[H1])-4.0/(1.0+exp(-(Y[V1]+27.0)/5.0))*Y[H1];   
F[N1] = 0.032*(Y[V1]+52.0)/(1.0-exp(-(Y[V1]+52.0)/5.0))*(1.0-Y[N1])-0.5*exp(-(Y[V1]+57.0)/40.0)*Y[N1];  
F[MN1] = (3.209*0.0001*(Y[V1]+30.0)/(1.0-exp(-(Y[V1]+30.0)/9.0)))*(1.0-Y[MN1]) + 3.209*0.0001*(Y[V1]+30.0)/(1.0-exp((Y[V1]+30.0)/9.0))*Y[MN1]; 
F[S1] = -Y[S1]/TAUSYN;  

current[I_NA2] = G_NA*Y[H2]*pow(Y[M2],3.0)*(Y[V2]-E_NA);
current[I_K2] = G_K*pow(Y[N2],4.0)*(Y[V2] - E_K);
current[I_M2] = G_M*Y[MN2]*(Y[V2] - E_K);
current[I_L2] =  G_L*(Y[V2] - E_L);
current[I_S2] =   gsyn2*Y[S1]*(Y[V2] - E_SYN) +G_SYN2*Y[S2]*(Y[V2]-E_SYN);
F[V2] = (I_APP + EPSILON - current[I_NA2] - current[I_K2] - current[I_M2] - current[I_L2] - current[I_S2])/CM;
F[M2] = 0.32*(Y[V2]+54.0)/(1.0-exp(-(Y[V2]+54.0)/4.0))*(1.0-Y[M2])-0.28*(Y[V2]+27.0)/(exp((Y[V2]+27.0)/5.0)-1.0)*Y[M2];   
F[H2] = 0.128*exp(-(Y[V2]+50.0)/18.0)*(1.0-Y[H2])-4.0/(1.0+exp(-(Y[V2]+27.0)/5.0))*Y[H2];   
F[N2] = 0.032*(Y[V2]+52.0)/(1.0-exp(-(Y[V2]+52.0)/5.0))*(1.0-Y[N2])-0.5*exp(-(Y[V2]+57.0)/40.0)*Y[N2];  
F[MN2] = (3.209*0.0001*(Y[V2]+30.0)/(1.0-exp(-(Y[V2]+30.0)/9.0)))*(1.0-Y[MN2]) + 3.209*0.0001*(Y[V2]+30.0)/(1.0-exp((Y[V2]+30.0)/9.0))*Y[MN2]; 
F[S2] = -Y[S2]/TAUSYN;  
return 0;
   }

void scan_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("state.data","r");
for(i=0;i<N;i++) fscanf(sp,"%lf\n",&Y[i]);
fclose(sp);}

void dump_(Y) 
double Y[N];
{FILE *fopen(),*sp;
int i;
sp = fopen("end.data","w");
for(i=0;i<N;i++) fprintf(sp,"%.16f\n",Y[i]);
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

