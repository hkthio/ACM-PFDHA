#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double interp(double *, double *, double, int);
float sort_down(int, double *, double *);
float sort_up(int, double *, double *);
int main()
   {
   double amp[100], rate[42], weight[42],cum[42],tot,fract,realweight[42];
   double percent[]={.05,.15,.5,.85,.95},mean;
   int idum,i,j,k,l,nEpist,nTarget;

   scanf("%d %d",&nEpist,&nTarget);

   /*Read in weights of th epistemic branches */
   tot=0.0;
   scanf("%d",&idum);
   for(j=0;j<nEpist;j++) 
      {
      scanf("%lf",(realweight+j));
      tot+=realweight[j];
      }
   for(j=0;j<nEpist;j++) realweight[j]/=tot;

   /*Loop over target amplitudes */
   printf("Displacement        Mean");
   for(k=0;k<5;k++) printf("    %5.3f    ",percent[k]);
   printf("\n");
   for(l=0;l<nTarget;l++)
      {
      mean=0.0;
      scanf("%lf",(amp+l));
      /*printf("Original\n");*/
      for(i=0;i<nEpist;i++)
         {
         scanf("%lf",(rate+i));
         weight[i]=realweight[i];
         mean+=weight[i]*rate[i];
         /*printf("%12.6e %12.6e\n",rate[i],weight[i]);*/
         }

      /*printf("Sorted:\n");*/
      sort_up(nEpist,rate,weight);
      /*for(i=0;i<nEpist;i++)printf("%12.6e %12.6e\n",rate[i],weight[i]);*/
      for(i=0;i<nEpist;i++) 
         {
         cum[i]=0.0;
         for(j=0;j<i;j++) cum[i]+=weight[j];
         }
      printf("%12.6e %12.6e ",amp[l],mean);
      for(k=0;k<5;k++)
         {
         fract=interp(cum,rate,percent[k],nEpist);
         printf(" %12.6e",fract);
         }
      printf("\n");
      }
   }






float sort_up(int n, double *a, double *w)
   {
   int i,j;
   double x,y;
   i=1;
   while(i<n)
      {
      j=i;
      while(j>0 && a[j-1] > a[j])
         {
         x=a[j];
         y=w[j];
         a[j]=a[j-1];
         w[j]=w[j-1];
         a[j-1]=x;
         w[j-1]=y;
         j-=1;
         }
      i++;
      }
   }

float sort_down(int n, double *a, double *w)
   {
   int i,j;
   float x,y;
   for(i=1;i<n;i++)
      {
      x=a[i];
      y=w[i];
      j=i-1;
      while(j >=0 && a[j] < x)
         {
         a[j+1]=a[j];
         w[j+1]=w[j];
         j--;
         }
      a[j+1]=x;
      w[j+1]=y;
      }
   }



double interp(double *x,double *y,double x0,int nx)
   {
/*    linear interpolation/extrapolation. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996 */
   int i,j,isign;
   double interp;
   j=0;
   isign=-1;
   if(x[0] <  x[nx-1]) isign=1;
   for(i=0;i<nx;i++) 
      {
      if(isign*(x0-x[i]) >  .0) j++;
      /*printf("%d %d %lf %lf %lf\n",i,j,x[i],y[i],isign*(x0-x[i]) );*/
      }
   /*printf("j=%5d\n",j);*/
   if(j >=  nx-1) 
      {
      interp=y[nx-1]+(y[nx-1]-y[nx-2])*(x0-x[nx-1])/(x[nx-1]-x[nx-2]);
      }
   else
      {
      if(j <  0) j=0;
      interp=y[j]+(y[j+1]-y[j])*(x0-x[j])/(x[j+1]-x[j]);
      }
   return interp;
   }

double interpl(double *xx,double *yy,double xx0, int nx)
   {
   /*    logarithmic interpolation of function y(x(i))
         where both x should be mononous and both x and y should
         be positive. datapoint outside the interval x(1)-x(nx)
         are extrapolated, using the nearest dy/dx.
         
                                     Hong Kie Thio, January 1996*/

   double x0,*x,*y,interpl;
   int i;

   x = (double *) calloc(nx, sizeof(double));
   y = (double *) calloc(nx, sizeof(double));

   for(i=0;i<nx;i++)
      {
      if(xx[i] <  .0  ||  yy[i] <   .0) 
         {
         printf("%f %f Negative x or y, cannot take log\n",xx[i],yy[i]);
         return -1;
         }
      else
         {
         if(yy[i] == 0.0) yy[i]=1.e-308;
         if(xx[i] == 0.0) xx[i]=1.e-308;
         x[i]=log10(xx[i]);
         y[i]=log10(yy[i]);
         }
      }
   x0=log10(xx0);
   interpl=interp(x,y,x0,nx);
   if(interpl <  -308.)
      {
      interpl=-308.;
      printf("Underflow, set to 1e-308\n");
      }
   else if(interpl >  36)
      {
      interpl=308.;
      printf("Overflow, set to 1e308\n");
      }
   interpl=pow(10.,interpl);
   return interpl;
   }


