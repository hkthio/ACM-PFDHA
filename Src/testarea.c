#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "PHA.h"
int main()
   {
   double p,R,A[]= {625, 2500, 10000, 22500, 40000};
   int i,j,nArea=5;
   for(i=0;i<100;i++)
      {
      R=i*5;
      printf("%5.1f ",R);
      for(j=0;j<nArea;j++)
         {
         p=probSurfRupArea(R,A[j]);
         printf("%12.6e ",p);
         }
      printf("\n");
      }
   }

double probSurfRupArea(float R,float A)
   /* probability of surface rupture from Petersen et al. (2011) */
   {
   int i;
   double probSurfRup;
   double area[] =  {625, 2500, 10000, 22500, 40000 };
   double a[] =     {-1.1470, -0.9000, -1.0114, -1.0934, -1.1538 };
   double b[] =     { 2.1046,  0.9866,  2.5572,  3.5526,  4.2342 };
   double sigma[] = { 1.2508,  1.1470,  1.0917,  1.0188,  1.0177 };
   double r[5][3] = { { .001, 100., 200. }, { .001, 100., 200. }, { .001, 100., 200. }, { .001, 150., 300. }, { .001, 200., 400. } };
   double p[5][3] = { { 74.541, 7.8690, 2.0108}, {87.162, 4.8206, 2.6177 }, {90.173, 18.523, 6.6354} , { 87.394, 19.592, 7.0477} , { 92.483, 18.975, 7.4709} };
   i=0;
   while (i<5 && area[i]< A) i++;
   if(R<r[i][2])
      {
      probSurfRup=interpl(r[i],p[i],R,3)/100.;
      }
   else
      {
      probSurfRup=exp(a[i]*log(R)+b[i]);
      }
   return probSurfRup;
   }

