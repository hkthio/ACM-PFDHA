#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "PHA.h"
#include "../Sphere.h"
#include "../CatUtil.h"
#define DEBUG 0

/* test Displacement Functions */

int main(int ac, char **av)
   {
   enum SLP{AVE, MAX};
   FILE *fp;
   int i,j,k,l,flag;
   double mag=7.5, Dav, Dmx, lnDav, lnDmx, D;
   double Target[1000],lnTarget,Sigma;
   double p, loverL,R, Ldim;
   double A,rate,pR,pS,pA,pSurfRup,ptarget,p1,pp,pD[6][1000];
   double D0, D1, D2, D3, D4, D5;
   double D05,D15,D50,D85,D95;
   int iAcc,iCmplx;

   mag=7.5;
   rate=1./140;
   ptarget=1./475;
   Dav=SlipScaling(mag,0,AVE,&Sigma);
   Dmx=SlipScaling(mag,0,MAX,&Sigma);
   lnDav=log(Dav);
   lnDmx=log(Dmx);
   printf("M=%5.6f Dav=%7.3f Dmx=%7.3f\n",mag,Dav,Dmx);

   for(i=0;i<=20;i++)
      {
      loverL=i*.025;
      for(j=0;j<1000;j++)
         {
         Target[j]=.002+j*j*.0002;
         lnTarget=log(Target[j]);
         pD[0][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,6,&D0);
         pD[1][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,7,&D1);
         pD[2][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,10,&D2);
         pD[3][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,11,&D3);
         pD[4][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,12,&D4);
         pD[5][j]=probFD(mag,loverL,lnDav,lnDmx,lnTarget,13,&D5);
         }
      printf("D %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",loverL,D0/Dav,D1/Dmx,D2/Dav,D3/Dmx,D4/Dav,D5/Dmx);
      for(j=0;j<6;j++)
         {
         D95=interp(pD[j],Target,0.95,1000);
         D85=interp(pD[j],Target,0.85,1000);
         D50=interp(pD[j],Target,0.50,1000);
         D15=interp(pD[j],Target,0.15,1000);
         D05=interp(pD[j],Target,0.05,1000);
         printf("F%d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",j,loverL,D05,D15,D50,D85,D95);
         }
      }
   }
