/*********************************************************************72
      program pnEMG
c
c    This program implements the cumulative distribution function
c     P_nEMG(x | mu, sigma, nu) of the nEMG distribution.
c
c    1. Calculation is based on the property
c
c     P_nEMG(x | mu, sigma, nu) = 1 - P_EMG(-x | -mu, sigma, nu)
c
c     , where P_EMG is the CDF of EMG distribution.
c
c    2. P_EMG is implemented in subroutine 'p_EMG',
c       which is a Fortran translation of the R function gamlss.dist::pexGAUS.
c
c    3. The complement of cumulative distribution function is
c       1 - P_nEMG(x | mu, sigma, nu) = P_EMG(-x | -mu, sigma, nu)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Created:
c
c    May 17, 2022
c
c  Author:
c
c    Brian Chiou
c
c   Compile: 
c           gfortran pnEMG.f -ffixed-line-length-none -static -o pnEMG
c*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double pEMG(double, double, double, double);

int main(int ac,char **av)
   {
   int i;
   double x, mu, sigma, nu;
   double fM,fl2L,p;
   double c0, c1, cv1, cv2, cv3, cv5, cv6, cn, m1, m2, m3;
   double l2L, l2Lf, Mag, x_star, sigmaeq, sigmatot,y;

   c0=1.30182;
   m1=3.50698;
   m2=0.76397;
   m3=7.10;
   c1=1.37534;
   cv1=1.06811;
   cv2=-0.75483;
   cv3=0.24120;
   cv5=1.22339;
   cv6=-0.95650;
   cn=10.;

   Mag=8.0;
   l2L=0.05;

   if(l2L <= 0.5)
      {
      l2Lf=l2L;
      }
   else
      {
      l2Lf=1.-l2L;
      }
   fM=m2*(Mag-m3)+((m2-m1)/cn)*log((1.+exp(-cn*(Mag-m3)))/2.);
   x_star=sqrt(1.-pow(l2L-0.5,2)/pow(0.5,2));
   fl2L=c1*(x_star-1.);
   sigma=cv3;
   sigmaeq=fmax(cv1*exp(cv2*fmax(Mag-6.1,0.0)),0.4);
   sigmatot=sqrt(pow(sigma,2)+pow(sigmaeq,2));
   nu=cv5*exp(cv6*l2Lf);

   mu=-(c0+fM+fl2L);

   printf("             x         x_star             mu         sigmaeq          sigma             nu     P[X <= x]\n");
   for(i=1;i<100;i++)
      {
      x=-log(pow(10,i/10.-5.));
      p = 1 - pEMG(x, mu,sigmatot,nu);
      printf("%14.9lf %14.9lf %14.9lf %14.9lf %14.9lf %14.9lf %14.9lf \n",exp(-x),x_star, -mu, sigmaeq,sigmatot, nu, p);
      }
   }
