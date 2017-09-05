/* emulator for homogenous all to all network of N neurons*/
#include <stdio.h>
#include <math.h>
#define N 2
#define INCREMENT 100
#define SPIKES 50000
#define F2OFF 0
#define ELAPSED 0
#define INTERVAL 1
#define DEBUG 0

main()
{
FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;
fp1 = fopen("f11.data","r");
fp2 = fopen("nperiod.data","r");
fp3 = fopen("nphi.data","r");
fp4 = fopen("sumf21.data","r");
fp5 = fopen("f21.data","r");
fp6 = fopen("endphi.data","w");
fp7 = fopen("endf21.data","w");
int mark,i,k,j,n,min,count,first;
int flag, update[N];
double rem,dummy;
double f1[INCREMENT+1];
double f2[INCREMENT+1];
double period[N];
double time_to_next_spike[N];
double saved_time[N];
double phi[N];
double saved_phi[N];
double phi_mirror[N];
double sumf21[N];
double elapsed_time=0.0;
for (n=0;n<INCREMENT;n++)
  {
  fscanf(fp1,"%lf %lf\n",&dummy,&f1[n]);
  fscanf(fp5,"%lf %lf\n",&dummy,&f2[n]);
  if (F2OFF) f2[n]=0.0;
  }
 f1[INCREMENT]= 0.0; /* required by frame of reference)*/
/* if (f1[INCREMENT] < 0) f1[INCREMENT]=0.0;   causality */
f2[INCREMENT]=f1[0]; 
  if (F2OFF) f2[INCREMENT]=0.0;
f2[0]= 0.0;
for (n=0;n<N;n++)
  {
  fscanf(fp2,"%lf\n",&period[n]);
  fscanf(fp3,"%lf\n",&phi[n]);
  fscanf(fp4,"%lf\n",&sumf21[n]);
  }


for (n=0;n<N;n++) saved_time[n] = -1.0;

for(k=0;k<SPIKES;k++)
   {  
    if(DEBUG) for (n=0;n<N;n++) printf(" A %f \n",phi[n]);
     time_to_next_spike[0] = period[0]*(1.0 - phi[0]);
     min= 0;
     
/* Figure out which neuron (s) will spike next without further inputs*/

    for (n=1;n<N;n++)
     {
     time_to_next_spike[n] = period[n]*(1.0 - phi[n]);
     if(time_to_next_spike[n] < time_to_next_spike[min] ) min = n;
     }
        /*printf("%f %d\n",time_to_next_spike[min],min);*/
fflush(stdout);
    for (n=0;n<N;n++) 
         {update[n]=0;
         }
 
    count = 0;
        first = 1; 
/* count how many neurons will spike at that time and mark them in the update array */

    for (n=0;n<N;n++)
     {
     if(time_to_next_spike[n] - time_to_next_spike[min] <= 0.0)
       {update[n] = 1;
        if(INTERVAL&&first) printf("%d %f\n",n, time_to_next_spike[n]);
        if(!first&&INTERVAL) printf("%d %f\n",n,0.000000);
        first = 0;
          count = count + 1;}

     }
     elapsed_time = elapsed_time + time_to_next_spike[min];

/* store the f2 contributions of the next firing to all neurons except the one firing*/

    if(DEBUG) for (n=0;n<N;n++) printf(" B %f \n",phi[n]);
    for (n=0;n<N;n++)
     {
     if(update[n]==1) 
         {if(ELAPSED) printf("%f %d\n",elapsed_time,n);
          fflush(stdout);
           for(j=0;j<N;j++)
              {
              if (j!=n) 
               {phi_mirror[j] = phi[j] + time_to_next_spike[min]/period[j];
if(update[j]==1) phi_mirror[j] = 0.0;
                mark = (phi_mirror[j]*INCREMENT)/1;
                rem = (phi_mirror[j]*INCREMENT) - mark;
if(DEBUG) printf("3  k= %d n = %d j = %d mark = %d phi_mirror[j] = %.16f sumf21[j] = %f rem = %f \n", k,n,j,mark,phi_mirror[j],sumf21[j],rem); 
if( mark<0 || mark == INCREMENT) {mark = 0;
                           rem= 0.0;}
                sumf21[j] = sumf21[j] + f2[mark]  + rem*(f2[mark+1] -f2[mark]);
if(DEBUG) printf("4  k= %d n = %d j = %d mark = %d phi_mirror[j] = %f sumf21[j] = %f rem = %f \n", k,n,j,mark,phi_mirror[j],sumf21[j],rem); 
               }
              }

         }
    }

    if(DEBUG) for (n=0;n<N;n++) printf(" C %f \n",phi[n]);
/*update the f2 contributions from the previous cycle for the neuron (s) that is firing, plus f1 if other neurons are also firing */
/*if an input is received by a neuron right after it has just received one, don't use the current phase but rather the phase at the previous input */
    for (n=0;n<N;n++)
     {
     if(update[n]==1) 
         {phi[n]=0.0 - sumf21[n] - (count-1)*f1[0];
        if(count>1) saved_phi[n] =phi[n];
        if(count>1) saved_time[n] = elapsed_time;
if(DEBUG) { printf("1  %f %f ",phi[n],sumf21[n]);
         printf("%d %d %d",count,k,n);
          printf("\n");}
/*         if(phi[n] >1.0) phi[n]=1.0;*/
          sumf21[n]=0;
         }
      }
fflush(stdout);
    if(DEBUG) for (n=0;n<N;n++) printf(" D %f \n",phi[n]);

/* update the phase and add the f1 contributions for the neurons that aren't spiking*/

    for (n=0;n<N;n++)
     {
     if(update[n]==0) 
       {
flag = 0;
if(DEBUG) printf("7  k= %d n = %d phi[n] = %.16f  \n", k,n,phi[n]); 
         phi[n] = phi[n] + time_to_next_spike[min]/period[n];
          mark = (phi[n]*INCREMENT)/1;
          rem = (phi[n]*INCREMENT) - mark;
if(DEBUG) printf("6  k= %d n = %d mark = %d phi[n] = %.16f sumf21[n] = %f rem = %f \n", k,n,mark,phi[n],sumf21[n],rem); 
          if(mark==INCREMENT) phi[n] = phi[n] - count*f1[0];
          if(mark<INCREMENT && mark>=0) phi[n] = phi[n] - count*(f1[mark] + rem*(f1[mark+1] -f1[mark]));
          if(mark<0) phi[n] = phi[n] - count*f1[0] ; 
if(DEBUG) printf("5  k= %d n = %d mark = %d phi[n] = %.16f f1[mark] = %f rem = %f \n", k,n,mark,phi[n],f1[mark],rem); 
if(DEBUG) { printf("2  %f %f ",phi[n],f1[mark]);
         printf("%d %d %d",count,k,n);
          printf("\n");}
        }
    }
   }

for (n=0; n<N; n++) fprintf(fp6,"%.16f ",phi[n]);
for (n=0; n<N; n++) fprintf(fp7,"%.16f ",sumf21[n]);
fclose(fp1);
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp5);
fclose(fp6);
fclose(fp7);

}
