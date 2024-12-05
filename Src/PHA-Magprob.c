#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "PHA.h"
#include "../Sphere.h"
#include "../CatUtil.h"
#define DEBUG 0

/* module for computing displacement hazard */
/* this part of the PHA framework for probabilistic seismic, tsunami and fault dispalcement hazard analysis */
/* version 0.1 - initial bare-bones code with the main rupture relations */

double probMagCH(float, float, float, float, float, float, double, double );
double probMagMX(double, double, double, double, double, double);
double probMagGR(double, double, double, double, double);

/* test program */
int main(int ac, char **av)
   {
   int k,nMag;
   double Mag,cMag,Mag0,minMag,maxMag,dMag,b, sigma;
   double p,pmax,pchar, cum,cummax,cumchar;

   b=1.0;
   sigma=.25;
   minMag=1.0;
   cMag=6.5;
   maxMag=8.0;
   dMag=0.1;
   Mag0=minMag+dMag/2.;
   nMag=(9.5-Mag0)/dMag+1;
   cum=cummax=0.0;
   for(k=0;k<nMag;k++)
      {
      Mag=Mag0+k*dMag;
      p=probMagGR(Mag,minMag,maxMag,dMag, b);
      pmax=probMagMX(Mag,cMag,minMag,maxMag,dMag, sigma);
      cum+=p;
      cummax+=pmax;
      printf("%6.3f %12.6e %12.6e\n",Mag,p,pmax);
      }
   printf("0.0 %12.6e %12.6f\n",cum,cummax);
   }

double probMagGR(double Mag, double minMag, double maxMag, double dMag, double b)
   {
   /* magnitude probability for the GR model */
   double p;
   double Mag0, Mag1;
   double a;

   Mag0=Mag-0.5*dMag;
   Mag1=Mag+0.5*dMag;
   if(Mag0 > maxMag)
      {
      p=0.0;
      }
   else
      {
      a = 1./(1. - exp(-b*(maxMag - minMag)) );
      p = a*(exp(-b*(Mag0 - minMag)) - exp(-b*(Mag1 - minMag)) );
      }
   return p;
   }

double probMagMX(double Mag, double cMag, double minMag, double maxMag, double dMag, double sigma)
   {
   /* magnitude probability for maximum magnitude model */
   double p;
   float Mag0, Mag1;
   double L0, U0, L1, U1;
   Mag0=Mag-0.5*dMag;
   Mag1=Mag+0.5*dMag;
   if(Mag0 > maxMag)
      {
      p=0.0;
      }
   else
      {
      L0=(Mag0-cMag)/sigma;
      U0=(Mag1-cMag)/sigma;
      L1=(minMag-cMag)/sigma;
      U1=(maxMag-cMag)/sigma;
      p=(normCDF(L0)-normCDF(U0))/(normCDF(L1)-normCDF(U1));
      }
   return p;
   }

double probMagCH(float Mag, float minMagCH, float MagGR, float minMag, float maxMag, float dMag, double N, double b)
   /* magnitude probability for characteristic model */
   /* doesn't work yet */
   {
   double p;
   float Mag0, Mag1, mU, magU;

   Mag0=Mag-0.5*dMag;
   Mag1=Mag+0.5*dMag;

   if(Mag1 <= minMagCH)
      /* GR section */
      {
      p = exp(-b*(Mag0 - minMag)) - exp(-b *(Mag1 - minMag));
      }
   else if(Mag0 > minMagCH)
      /* char mag section */
      {
      if ( mU < magU )
         {
         p = b*(exp(-b*(MagGR-minMag))) * dMag;
         }
      else
         {
         p = b*(exp(-b*(MagGR-minMag))) * (maxMag-Mag0);
         }
      }
   else
      /* mixed section */
      {
      }
   return p;
   }
