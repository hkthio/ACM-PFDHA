#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "PHA.h"
#include "../Sphere.h"
#include "../CatUtil.h"
#define DEBUG 0

#define maxMagBin 16
#define maxEpsBin 8
#define maxTrg 100

/*---------------------------------------------------------------------------*/
/* module for computing displacement hazard                                  */
/* this part of the PHA framework for probabilistic seismic, tsunami and     */
/* fault displacement hazard analysis                                        */
/* version 0.1 - initial bare-bones code with the main rupture relations     */
/* version 0.2 - initial bare-bones code with the main and secondary rupture */
/*               relations and fractiles                                     */
/*                                                                           */
/* epistemic:                                                                */
/*   - probability of surface rupture                                        */
/*   - slip rate                                                             */
/* aleatory:                                                                 */
/*   - magnitude                                                             */
/*   - displacement                                                          */
/*   -                                                                       */
/*   -                                                                       */
/*   -                                                                       */
/*---------------------------------------------------------------------------*/

int getBin(int, double *, double );

int main(int ac, char **av)
   {
   enum SLP{AVE, MAX};
   FILE *fp, *fp1, *fp2;
   int i,j,k,l,flag,fdm;
   int hw=0,nmag,nrate[12],nsr,isr[11];
   double Wt,mag[12],magWt[12], Dav, Dmx, lnDav, lnDmx, D,DD,logD,sigD;
   double lnTarget,Sigma,srWt[11];
   double p1[maxTrg][6][maxTrg][5],p2[maxTrg][6][maxTrg][5],pp[maxTrg][5][maxTrg][6];
   double p3[maxTrg],p4[maxTrg],p5[maxTrg];
   double pk[12][11][maxTrg], pm[12][maxTrg], pl[12][maxTrg];
   double loverL,R, Ldim;
   double A,rate[12][3],rateWt[12][3];
   double pR_main,pR_sec,pS,pA_main,pA_sec,pSurfRup,ptarget, pD_main,pD_sec;
   double sliprate,fdmWt;
   double epsilon;
   /* bins */
   double epsBin[]={-3.,-2.,-1.,0.,1.,2.,3};
   double magBin[]={6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5,8.75,9.0,9.25};
   double EpsMagBin[maxTrg][maxEpsBin][maxMagBin]={0.0};
   int nEpsBin=7,nMagBin=14,ieps,imag;
   double retp[] = {72, 475, 975, 2475, 5000, 10000, 25000, 50000, 100000, 250000, 
             500000, 1000000, 2500000, 5000000, 10000000};
   double Target[] = { 1.000E-03,1.205E-03,1.451E-03,1.592E-03,1.748E-03,1.918E-03,
             2.105E-03,2.536E-03,2.783E-03,3.055E-03,3.353E-03,3.680E-03,4.038E-03,
             4.432E-03,4.864E-03,5.339E-03,5.860E-03,6.431E-03,7.058E-03,7.747E-03,
             8.502E-03,9.331E-03,1.024E-02,1.124E-02,1.234E-02,1.354E-02,1.486E-02,
             1.631E-02,1.790E-02,1.964E-02,2.156E-02,2.366E-02,2.597E-02,2.850E-02,
             3.128E-02,3.433E-02,3.768E-02,4.136E-02,4.539E-02,4.982E-02,5.468E-02,
             6.001E-02,6.586E-02,7.228E-02,7.933E-02,8.707E-02,9.556E-02,1.049E-01,
             1.151E-01,1.263E-01,1.387E-01,1.522E-01,1.670E-01,1.833E-01,2.012E-01,
             2.208E-01,2.423E-01,2.660E-01,2.919E-01,3.204E-01,3.516E-01,3.859E-01,
             4.236E-01,4.649E-01,5.102E-01,5.599E-01,6.146E-01,6.745E-01,7.403E-01,
             8.125E-01,8.917E-01,9.787E-01,1.074E+00,1.179E+00,1.294E+00,1.420E+00,
             1.558E+00,1.710E+00,1.877E+00,2.060E+00,2.261E+00,2.482E+00,2.724E+00,
             2.989E+00,3.281E+00,3.601E+00,3.952E+00,4.338E+00,4.761E+00,5.225E+00,
             5.734E+00,6.294E+00,6.908E+00,7.581E+00,8.321E+00,9.132E+00,1.002E+01,
             1.200E+01,1.500E+01,2.000E+01 };
   int iAcc,iCmplx, nret=1 ;

   /* read geometry */
   scanf("%d %lf",&fdm,&fdmWt);
   scanf("%lf",&loverL);
   scanf("%lf",&Ldim);
   scanf("%lf",&R);
   scanf("%d",&iAcc);
   scanf("%d",&iCmplx);
   scanf("%d",&nsr);
   /* read probability of surface rupture branches */
   for(k=0;k<nsr;k++) 
      {
      scanf("%d %lf",(isr+k),(srWt+k));
      }

   /* read magnitude and rate branches */
   scanf("%d",&nmag);
   for(i=0;i<nmag;i++) 
      {
      scanf("%lf %lf",(mag+i),(magWt+i));
      scanf("%d",(nrate+i));
      for(j=0;j<nrate[i];j++) scanf("%lf %lf",(rate[i]+j),(rateWt[i]+j));
      }
   A=Ldim*Ldim;

   printf("%6.4f %7.2f %8.2f %d %d\n",loverL,Ldim,R,iAcc,iCmplx);
   /* read magnitude and rate branches */
   printf("%d\n",nmag);
   for(i=0;i<nmag;i++) 
      {
      printf("%2d Mag: %5.3f Wt: %6.4f\n",i,mag[i],magWt[i]);
      for(j=0;j<nrate[i];j++) 
                        printf("Rates %2d %12.6e %6.4f\n",j,rate[i][j],rateWt[i][j]);
      }
   /* end input */

   for(j=0;j<100;j++)
      {
      p3[j]=0.0;
      p4[j]=0.0;
      p5[j]=0.0;
      Target[j]=3.*Target[j];
      }
   for(k=0;k< nsr;k++)
      {
      for(i=0;i<nmag;i++)
         {
         pS=probSurfRup( mag[i],isr[k]);
         Dav=SlipScaling(mag[i],3,AVE,&Sigma);
         Dmx=SlipScaling(mag[i],3,MAX,&Sigma);
         printf("M= %4.2f Dav = %6.2f, Dmax = %6.2f\n",mag[i],Dav,Dmx);
         for(l=0;l<nrate[i];l++)
            {
            if(rate[i][l] > 0)
               {
               printf("Return period = %10.2f years\n",1./rate[i][l]);
               }
            else
               {
               printf("Return period = %10.2f years\n",Dav/(-.001*rate[i][l]));
               }
            }
         lnDav=log(Dav);
         lnDmx=log(Dmx);
         pA_main=probSurfRupArea(0,A,mag[i],fdm,hw);
         pA_sec=probSurfRupArea(R,A,mag[i],fdm,hw);
         pR_main = probSurfRupR(R,Ldim, iAcc, iCmplx, fdm);
         pR_sec  = probSurfRupR(R,Ldim, iAcc, iCmplx, fdm);
         pSurfRup = pS*pR_sec*pA_sec;
         printf("  iAcc,         R,        pS,   pR_main,    pR_sec,   pA_main,    pA_sec,         D,         d\n");
         printf("%6d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e ",
                iAcc,R,pS,pR_main,pR_sec,pA_main,pA_sec);
         flag=0;
         lnTarget=1.;
         pD_main=probFD(mag[i],loverL,lnDav,lnDmx,lnTarget,fdm,&logD, &sigD);
         printf("%10.4e ",pD_main);
         pD_sec=probFDDistr(mag[i],R,lnDav,lnDmx,lnTarget,fdm,&D,hw);
         printf("%10.4e",pD_sec);
         printf("\n");
         
         for(j=0;j<100;j++)
            {
            /*Target[j]=pow(10,j*.015-15);*/
            lnTarget=log(Target[j]);
            pD_main = probFD(mag[i],loverL,lnDav,lnDmx,lnTarget,fdm,&logD,&sigD);
            epsilon=(lnTarget-logD)/sigD;
            pD_sec  = probFDDistr(mag[i],R,lnDav,lnDmx,lnTarget,fdm,&D,hw);
            for(l=0;l<nrate[i];l++)
               {
               Wt=srWt[k]*rateWt[i][l]*magWt[i];
               if(rate[i][l] > 0) 
                  {
                  p1[i][l][j][k] = rate[i][l]*pS*pD_main*pR_main*pA_main;
                  p2[i][l][j][k] = rate[i][l]*pS*pD_sec*pA_sec;
                  pp[i][l][j][k] = rate[i][l]*pS*pD_main;
                  }
               else
                  {
                  p1[i][l][j][k] = (-.001*rate[i][l]/Dav)*pS*pD_main*pR_main*pA_main;
                  p2[i][l][j][k] = (-.001*rate[i][l]/Dav)*pS*pD_sec*pA_sec;
                  pp[i][l][j][k] = (-.001*rate[i][l]/Dav)*pS*pD_main;
                  }
               p3[j] += Wt*p1[i][l][j][k];
               p4[j] += Wt*p2[i][l][j][k];
               p5[j] += Wt*pp[i][l][j][k];

               /* fill disaggregation bins */
               ieps=getBin(nEpsBin,epsBin,epsilon);
               imag=getBin(nMagBin,magBin,mag[i]);
               printf("Epsilon = %f, ieps = %d, Mag = %f, imag = %d\n",
                      epsilon, ieps, mag[i], imag);
               EpsMagBin[j][ieps][imag]+=Wt*pp[i][l][j][k];
               }
            }
       /*for(j=0;j<nret;j++)
            {
            D=interpl(p1,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            D=interpl(p2,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            D=interpl(p3,Target,1./retp[j],1000);
            printf("%10.4e,",D);
            }
         printf("\n");*/
         } /* end Mag loop */
      }    /* end pSurfRup loop */
   fp=fopen("o_PFDHA","w");
   for(j=0;j<100;j++)
      {
      fprintf(fp,"%12.6e ",Target[j]);
      /*
      for(i=0;i<nmag;i++)
         {
         for(l=0;l<nrate[i];l++) printf("%12.6e,%12.6e,",p1[i][l][j],p2[i][l][j]);
         }*/
      fprintf(fp,"%12.6e %12.6e %12.6e %9.3f %12.6e\n",
              p3[j],p4[j],p3[j]+p4[j],1./(p3[j]+p4[j]),p5[j]);
      }
   fclose(fp);


   fp=fopen("o_deagg.csv","w");
   fprintf(fp,"%2d,%2d\n",nsr,nmag);
   fprintf(fp,"pR, prWt, Mag , magWt,  Rate     ,rateWt,");
   for(j=0;j<100;j++) fprintf(fp,"%12.6e,",Target[j]);
   fprintf(fp,"\n ");
   /*first loop over the epistemic branches to compute the hazard curves*/
   for(k=0;k< nsr;k++)
      {
      for(i=0;i<nmag;i++)
         {
         for(l=0;l<nrate[i];l++)
            {
            fprintf(fp,"%2d,%12.6e,%5.3f,%12.6e,%12.6e,%12.6e,",
                    isr[k],srWt[k],mag[i],magWt[i],rate[i][l],rateWt[i][l]);
            for(j=0;j<100;j++)
               {
               fprintf(fp,"%12.6e,",p1[i][l][j][k]);
               /* epistemic hazard curves */
               pk[k][l][j] += magWt[i]*p1[i][l][j][k];
               /* magnitude branch (non-epistemic for now) */
               pm[i][j] += p1[i][l][j][k];
               /* rate branch */
               pl[l][j] += p1[i][l][j][k];
               }
            fprintf(fp,"\n");
            }
         }
      }
   fclose(fp);

   /* write out marginal hazard per probSurfrup function */
   fp =fopen("o_epistemic.xyz","w");
   fp2=fopen("o_disag-eps-mag.csv","w");
   fprintf(fp,"%3d\n",maxTrg);
   fprintf(fp,"        %4d",nsr*nrate[0]);
   for(k=0;k< nsr;k++)      
      {
      for(l=0;l< nrate[0];l++) fprintf(fp," %12.6e",fdmWt*srWt[k]*rateWt[0][l]);
      }
   fprintf(fp,"\n");
   for(j=0;j<100;j++) 
      {
      fprintf(fp, "%12.6e",Target[j]);
      for(k=0;k< nsr;k++)      
         {
         for(l=0;l< nrate[0];l++) fprintf(fp," %12.6e",pk[k][l][j]);
         }
      for(i=0;i < nMagBin;i++) 
         {
         fprintf(fp2,"%7.4f,%5.3f",Target[j],magBin[i]);
         for(l=0;l <=nEpsBin;l++) fprintf(fp2,",%12.6e",EpsMagBin[j][l][i]);
         fprintf(fp2,",%12.6e\n",p5[j]);
         }
      fprintf(fp,"\n");
      }
   fclose(fp);
   fclose(fp1);
   fclose(fp2);
   }


int getBin(int nbins, double *bins, double value)
   {
   int i=0;
   while(i<nbins && value>bins[i]) i++;
   return i;
   }
