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
   int hang;
   double mag=7.5, Dav, Dmx, lnDav, lnDmx, D,DD;
   double Target[1000],lnTarget,Sigma;
   double p1[1000],p2[1000],p3[1000], loverL,R, Ldim;
   double A,rate,pR_main,pR_sec,pS,pA_main,pA_sec,pSurfRup,ptarget,pp;
   double pD_main,pD_sec;
   double sliprate;
   /*double retp[] = {2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000};*/
   double retp[] = {475, 975, 2475, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000};
   int iAcc,iCmplx, nret=1 ;

   mag=7.0;
   sliprate=.01;
   Dav=SlipScaling(mag,0,AVE,&Sigma);
   Dmx=SlipScaling(mag,0,MAX,&Sigma);
   lnDav=log(Dav);
   lnDmx=log(Dmx);
   rate=sliprate*1.e-3/Dav;
   rate=1./140.;

   Ldim=200.;
   A=Ldim*Ldim;
   iAcc=0;
   iCmplx=0;

   pS=probSurfRup( mag, 1);
   printf("Magnitude = %4.2f Sliprate = %6.2f Eventrate = %12.6e 1/yr Slip = %6.3f  pS=%8.6f Area=%9.2f Accuracy=%1d Complexity=%d\n",mag,sliprate,rate,Dav,pS,A,iAcc,iCmplx);
   printf("     1,");
   for(j=2;j<10;j ++) printf("  %8d,",j);
   for(j=0;j<nret;j++) printf("  %8d,  %8d,  %8d",10+3*j,11+3*j,12+3*j);
   printf("\n");
   printf("      ,          ,          ,          ,          ,          ,          ,");
   for(j=0;j<nret;j++) printf("  %8.0f,  %8.0f,",retp[j],retp[j]);
   printf("\n");
   printf("  iAcc,         R,        pS,   pR_main,    pR_sec,   pA_main,    pA_sec,         D,         d,");
   for(j=0;j<nret;j++) printf("    D-main,  D-second,     D-all");
   printf("\n");
   for(i=0;i<50;i++)
      {
      R=i*10;
      pA_main=probSurfRupArea(3,A);
      pA_sec=probSurfRupArea(R,A);
      k=6;
      for(l=1;l<5;l++)
         {
         loverL=0.4;
         iAcc=l;
         pR_main = probSurfRupR(R,Ldim, iAcc,iCmplx);
         pR_sec  = probSurfRupR(R,Ldim, iAcc, iCmplx);
         pSurfRup = pS*pR_sec*pA_sec;
         printf("%6d,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,%10.4e,",iAcc,R,pS,pR_main,pR_sec,pA_main,pA_sec);
         flag=0;
         lnTarget=1.;
         pD_main=probFD(mag,loverL,lnDav,lnDmx,lnTarget,k,&D);
         printf("%10.4e,",D);
         pD_sec=probFDDistr(mag,R,lnDav,lnTarget,k,&D,hang);
         printf("%10.4e,",D);
         for(j=0;j<1000;j++)
            {
            Target[j]=pow(10,j*.015-15);
            lnTarget=log(Target[j]);
            pD_main=probFD(mag,loverL,lnDav,lnDmx,lnTarget,k,&D);
            pD_sec=probFDDistr(mag,R,lnDav,lnTarget,k,&D,hang);
            p1[j]=rate*pD_main*pR_main*pA_main;
            p2[j]=rate*pD_sec*pA_sec;
            p3[j]=p1[j]+p2[j];
            /*printf("%2d %2d %6.3f %6.3f %6.3f %6.3f %6.3f %12.6e\n",i,k,loverL,Dav,Dmx,D,Target,p);*/
            /*printf("%2d %2d %6.3f %6.3f %6.3f %12.6e %12.6e\n",j,k,R,Dav,Target[j],pD_main,pD_sec);*/
            }
         for(j=0;j<nret;j++)
            {
            D=interpl(p1,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            D=interpl(p2,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            D=interpl(p3,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            }
         printf("\n");
         }
      }
   }
