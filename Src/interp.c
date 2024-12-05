#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int main(int ac, char **av)
   {
   double interpl(double *,double *,double,int);
   double interp(double *,double *,double,int);
   double a[1000],b[1000];
   double target;
   int i,ni;
   double result;

   sscanf(av[1],"%lf",&target);
   
   i=0;
   while(scanf("%lf %lf",(a+i),(b+i)) == 2) i++;
   ni=i;
   /*printf("%d pairs read\n",ni);*/

   result = interp(b,a,target,ni); 
   printf("%6.3lf\n",result);
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


