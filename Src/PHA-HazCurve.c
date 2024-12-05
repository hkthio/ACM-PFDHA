#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "PHA.h"
#include "../Sphere.h"
#include "../CatUtil.h"
#define DEBUG 0

/* module for computing displacement hazard */
/* this part of the PHA framework for probabilistic seismic, tsunami and fault displacement hazard analysis */
/* version 0.1 - initial bare-bones code with the main rupture relations */

/* test program */
int main(int ac, char **av)
   {
   enum SLP{AVE, MAX};
   FILE *fp;
   int i,j,k,l,flag;
   double mag=7.5, Dav, Dmx, lnDav, lnDmx, D;
   double Target[1000],lnTarget,Sigma;
   double p1[1000],p2[1000], loverL,R, Ldim;
   double A,rate,pR_main,pR_sec,pS,pA_main,pA_sec,pSurfRup,ptarget,pp;
   double pD_main,pD_sec;
   double sliprate;
   /*double retp[] = {2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000};*/
   double retp[] = {475, 975, 2475, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000};
   int iAcc,iCmplx, nret=9 ;

   mag=7.0;
   sliprate=10;
   Dav=SlipScaling(mag,0,AVE,&Sigma);
   Dmx=SlipScaling(mag,0,MAX,&Sigma);
   lnDav=log(Dav);
   lnDmx=log(Dmx);
   rate=sliprate*1.e-3/Dav;
   rate=1./140.;

   Ldim=25.;
   A=Ldim*Ldim;
   iAcc=1;
   iCmplx=0;

   pS=probSurfRup( mag, 1);
   printf("Magnitude = %4.2f Sliprate = %6.2f Eventrate = %12.6e 1/yr Slip = %6.3f  pS=%8.6f Area=%9.2f Accuracy=%1d Complexity=%d\n",mag,sliprate,rate,Dav,pS,A,iAcc,iCmplx);
   printf("     1,");
   for(j=2;j<8;j ++) printf("  %8d,",j);
   for(j=0;j<nret;j++) printf("  %8d,  %8d,",8+2*j,9+2*j);
   printf("\n");
   printf("      ,          ,          ,          ,          ,          ,          ,");
   for(j=0;j<nret;j++) printf("  %8.0f,  %8.0f,",retp[j],retp[j]);
   printf("\n");
   printf("loverL,         R,        pS,   pR_main,    pR_sec,   pA_main,    pA_sec,");
   for(j=0;j<nret;j++) printf("    D-main,  D-second,");
   printf("\n");
   for(i=0;i<1;i++)
      {
      R=pow(10,i*.1-2);
      R=64.;
      pA_main=probSurfRupArea(0,A);
      pA_sec=probSurfRupArea(R,A);
      k=0;
      for(l=0;l<1;l++)
         {
         loverL=l*.1;
         loverL=.5;
         pR_main = probSurfRupR(R,-1.,iAcc,iCmplx);
         pR_sec  = probSurfRupR(R,Ldim, iAcc, iCmplx);
         pSurfRup = pS*pR_sec*pA_sec;
         printf("%6.4f,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,",loverL,R,pS,pR_main,pR_sec,pA_main,pA_sec);
         flag=0;
         for(j=0;j<1000;j++)
            {
            Target[j]=pow(10,j*.015-14);
            lnTarget=log(Target[j]);
            pD_main=probFD(mag,loverL,lnDav,lnDmx,lnTarget,k,&D);
            pD_sec=probFDDistr(mag,R,lnDav,lnTarget,k,&D);
            p1[j]=rate*pD_main;
            p2[j]=rate*pD_sec;
            printf("%2d %2d %6.3f %6.3f %6.3f %10.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",i,k,loverL,Dav,Dmx,Target[j],rate,pA_main,pA_sec,pR_main,pR_sec,pD_main,pD_sec);
            }
         /*for(j=0;j<nret;j++)
            {
            D=interpl(p1,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            D=interpl(p2,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            }
         printf("\n");*/
         }
      }
   }
