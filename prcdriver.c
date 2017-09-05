#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lts.h"
#include "prclts.c"
#include "cblock.h"
#define SAMPLE 0


/* global variable declarations */

/* double F[N];
double Y[N];
long int *np;
double *xp; */
 
double gsyn,gsyn1=0.0,gsyn2=0.0,current[C],vset1,initial[N],period[2];
double ss1x,ss2x;
int flag1=0,flag2=0,flag3=0,flag4=0;
double f1,f2,f3;

void main(int argc, char *argv[])
{
double state[N];
double fstate[N];
double rtol[N],atol[N];

int itol,ijac,mljac,mujac;
int imas,mlmas,mumas;
int iout,lwork,liwork,lrcont;
int iwork[LIWORK];
double work[LWORK];
extern int sout_(),out_(),oldout_(),altout_(),altout2_();
extern int debugout_(),debugout2_();
extern int flag1,flag2,flag3,flag4;
extern double f1,f2,f3;
int n,idid,counter=0;

double x,xend,h;
FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*pi,*ip;
int i,j;
if(argc >1 ) gsyn = strtod(argv[1],NULL);
else gsyn = G_SYN;
fp1 = fopen("f11.data","w");
fp2 = fopen("f21.data","w");
fp3 = fopen("f31.data","w");
fp4 = fopen("f12.data","w");
fp5 = fopen("f22.data","w");
fp6 = fopen("f32.data","w");
pi = fopen("period.data","w");
ip = fopen("initial.data","w");
n=N;
for(i=0;i<7;i++) {iwork[i]=0;
                  work[i]=0.0;}
for(i=1;i<N;i++) {rtol[i]=1.0e-8; 
                  atol[i]=1.0e-12;}
rtol[0]=1.0e-8;
atol[0]=1.0e-10;
rtol[6]=1.0e-8;
atol[6]=1.0e-10;      

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
xend=ENDTIME;
gsyn1 = 0.0;
gsyn2 = 0.0;
radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        out_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
  gsyn1=gsyn;
for(counter=0;counter<INCREMENT;counter++)
 {x = 0.0;
  xend=  (counter*period[0])/(1.0*INCREMENT) ;
  state[0]=THRESHOLD;;
  state[1]=initial[1];
  state[2]=initial[2];
  state[3]=initial[3];
  state[4]=initial[4];
  state[11]=0.000000;
  if(counter>0) radau5_(&n,deriv1_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
  state[6]=THRESHOLD;
  state[7]=initial[7];
  state[8]=initial[8];
  state[9]=initial[9];
  state[10]=initial[10];
  state[5]=0.000000;
  h=1e-4;
  flag1=0;
  flag2=0;
  flag3=0;
  flag4=0;
  xend = xend + period[1];
  radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
  h=1e-4;
   do{
      xend = xend + 10.0;
      radau5_(&n,deriv1_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
        } while (!flag3);
  fprintf(fp1,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f1);
  fprintf(fp2,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f2);
  fprintf(fp3,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f3);
 /*printf("%f %f\n",period[0]*(f2 + (1.0*counter)/(1.0*INCREMENT)), period[0]*(1 + f1 - (1.0*counter)/(1.0*INCREMENT))); */
 }
  /*gsyn2=gsyn;
gsyn1 = 0.0;
for(counter=0;counter<INCREMENT;counter++)
 {x = 0.0;
  xend=  (counter*period[1])/(1.0*INCREMENT) ;
  state[6]=THRESHOLD;;
  state[7]=initial[7];
  state[8]=initial[8];
  state[9]=initial[9];
  state[10]=initial[10];
  state[5]=0.000000;
  if(counter>0) radau5_(&n,deriv1_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout2_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
  state[0]=THRESHOLD;
  state[1]=initial[1];
  state[2]=initial[2];
  state[3]=initial[3];
  state[4]=initial[4];
  state[11]=0.000000;
  h=1e-4;
  flag1=0;
  flag2=0;
  flag3=0;
  flag4=0;
  xend = xend + period[0];
  radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout2_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
  h=1e-4;
   do{
      xend = xend + 10.0;
      radau5_(&n,deriv1_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        altout2_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);
        } while (!flag3);
  fprintf(fp4,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f1);
fflush(fp4);
  fprintf(fp5,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f2);
  fprintf(fp6,"%f %f\n",(1.0*counter)/(1.0*INCREMENT),f3);
 }*/
/* This part samples the limit cycle evenly commented out for now*/
if(SAMPLE) {  xend = xend + 2*period[0];
  radau5_(&n,deriv_,&x,state,&xend,&h,
        rtol,atol,&itol,
        dummy,&ijac,&mljac,&mujac,
        mas,&imas,&mlmas,&mumas,
        sout_,&iout,
        work,&lwork,iwork,&liwork,&lrcont,&idid);}
  fprintf(pi,"%f %f\n",period[0],period[1]);
  for(n=0;n<N;n++) fprintf(ip,"%f\n",initial[n]);
fflush(stdout);
fclose(fp1);
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp5);
fclose(fp6);
fclose(pi);
fclose(ip);
}

