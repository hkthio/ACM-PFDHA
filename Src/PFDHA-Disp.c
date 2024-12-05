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

float SlipScaling(double mag, int i, int j, double *Sigma)
   {
   float Slip;
   /* Slip values from various empirical relations */
   /* Second index j */
   /* 0 Ave Slip */
   /* 1 Max Slip */
   /* First index i */
   /* 0 Wells & Coppersmith - SS */ 
   /* 1 Wells & Coppersmith - RV, poor*/
   /* 2 Wells & Coppersmith - N  */
   /* 3 Wells & Coppersmith - All*/
   int nscal = 4;
   float a[][2] = {-6.32,-7.03,-0.74,-1.84,-4.45,-5.90,-4.80,-5.46 };
   float b[][2] = { 0.90, 1.03, 0.08, 0.29, 0.63, 0.89, 0.69, 0.82 };
   float s[][2] = { 0.28, 0.34, 0.38, 0.42, 0.33, 0.38, 0.36, 0.42 };
   /*printf("%d %d %5.2f %5.2f %5.2f\n",i,j,a[i][j],b[i][j],s[i][j]);*/
   if(i >= nscal)
      {
      printf("No slip scaling relations: %d %d\n",i,nscal);
      exit(-1);
      }
   else
      {
      Slip=pow(10.,a[i][j]+b[i][j]*mag);
      *Sigma=s[i][j];
      return Slip;
      }
   return Slip;
   }

/* Probability of surface rupture */
struct FDEq *InitEq()
   {
   /* initialize the displacement equations */
   /* 1-6 = Petersen */
   int i,neq = 15;
   struct FDEq *D = malloc(sizeof(struct FDEq) * neq);
   for(i=0;i<6;i++)      
      {
      D[i].type=NORM;
      D[i].maperr=GAUSS;
      /*D[i].maperr=DBLEXP;*/
      }
   for(i=6;i< 14;i+=2)   D[i].type=GAMMA;
   for(i=7;i< 14;i+=2)   D[i].type=BETA;
   D[14].type=PEMG;
   D[14].maperr=GAUSS;
   /*D[14].maperr=DBLEXP;*/

   D[0].subtype=BIL;     D[1].subtype=QUAD;     D[2].subtype=ELLI;
   D[3].subtype=BILN;    D[4].subtype=QUADN;    D[5].subtype=ELLIN;
   D[6].subtype=TAKAO;   D[7].subtype=TAKAO;
   D[8].subtype=JNES;    D[9].subtype=JNES;
   D[10].subtype=MOSS;   D[11].subtype=MOSS;
   D[12].subtype=YOUNGS; D[13].subtype=YOUNGS;
   D[14].subtype=CHIOU;
   return D;
   }

double probFD(float mag, double loverL, double logD_av, double logD_max, double logTarget, int FDM, double *logD, double *sigD)
   /* function to compute the probabilistic fault displacement */
   /* FDM is the empirical relation */
   {
   enum EquType type;
   enum EquSubType subtype;
   float a,a1,a2,a3,a4,b,b1,b2,b3,b4,c1,c2,c3;
   float loverLp,sigTrunc=12.0;
   double x;
   double *Disp,probFD;
   struct FDEq *D;

   /* Chiou et al exra variables */
   double  mu, sigma, nu;
   double fM,fl2L;
   double c0, cv1, cv2, cv3, cv5, cv6, cn, m1, m2, m3;
   double  l2Lf, x_star, sigmaeq, sigmatot;



   D=InitEq();

   type=D[FDM].type;
   subtype=D[FDM].subtype;
   switch(type)
      {
      case NORM:
         if(DEBUG) printf("Normal distribution\n");
         switch(subtype)
            {
            case BIL:
               if(DEBUG) printf("Bilinear\n");
               a1=1.7969; a2=1.7658; b=8.5206; c1=-10.2855; c2=-7.8962;
               loverLp=((a2-a1)*mag+(c2-c1))/b;
               if(loverL < loverLp)
                  {
                  *logD=a1*mag+b*loverL+c1-log(100.);
                  *sigD=1.2906;
                  }
               else
                  {
                  *logD=a2*mag+c2-log(100.);
                  *sigD=0.9624;
                  }
               break;
            case QUAD:
               if(DEBUG) printf("Quadratic\n");
               *logD=1.7895*mag+14.4696*loverL-20.1723*pow(loverL,2)-10.54512-log(100.);
               *sigD=1.1346;
               break;
            case ELLI:
               if(DEBUG) printf("Elliptic\n");
               *logD=3.3041*sqrt(1.-4*pow((loverL-0.5),2))+1.7927*mag-11.2912-log(100.);
               *sigD=1.1348;
               break;
            case BILN:
               if(DEBUG) printf("Bilinear normalized\n");
               if(loverL < 0.3008)
                  {
                  *logD=8.2525*loverL-2.3010+logD_av;
                  *sigD=1.2962;
                  }
               else
                  {
                  *logD=0.1816+logD_av;
                  *sigD=1.0013;
                  }
               break;
            case QUADN:
               if(DEBUG) printf("Quadratic normalized\n");
               *logD=14.2824*loverL-19.8833*pow(loverL,2)-2.6279+logD_av;
               *sigD=1.1419;
               break;
            case ELLIN:
               if(DEBUG) printf("Elliptic normalized\n");
               *logD=3.2699*sqrt(1.-4*pow((loverL-0.5),2))-3.2749+logD_av;
               *sigD=1.1419;
               break;
            default: 
               printf("FD sub relation type not know: %d\n",subtype);
               return -1.;
            }
         probFD = pXceed ( *logD, logTarget,*sigD, sigTrunc);
         return probFD;
         break;
      case PEMG:
         switch(subtype)
            {
            case CHIOU:
               if(DEBUG) printf("Chiou et al. 2023\n");
               c0=1.30182; m1=3.50698; m2=0.76397; m3=7.10; c1=1.37534; cv1=1.06811; cv2=-0.75483; cv3=0.24120; cv5=1.22339; cv6=-0.95650; cn=10.;
               if(loverL <= 0.5)
                  {
                  l2Lf=loverL;
                  }
               else
                  {
                  l2Lf=1.-loverL;
                  }
               fM=m2*(mag-m3)+((m2-m1)/cn)*log((1.+exp(-cn*(mag-m3)))/2.);
               x_star=sqrt(1.-pow(loverL-0.5,2)/pow(0.5,2));
               fl2L=c1*(x_star-1.);
               sigma=cv3;
               sigmaeq=fmax(cv1*exp(cv2*fmax(mag-6.1,0.0)),0.4);

               sigmatot=sqrt(pow(sigma,2)+pow(sigmaeq,2));
               nu=cv5*exp(cv6*l2Lf);
               mu=(c0+fM+fl2L);
               *sigD=sigmatot;
               *logD=mu-nu;
               break;
            default: 
               printf("FD sub relation type not know: %d\n",subtype);
               return -1.;
            }
         /*printf("\nRunning pEMG: logTarget = %f, mu = %f, sigmatot = %f, nu = %f\n",logTarget, mu, sigmatot, nu);*/
         probFD = pEMG(-logTarget, -mu,sigmatot,nu);
         /*printf("probFD = %f\n",probFD);*/
         if(probFD > 1.0 ) printf("Error pEMG = %f\n",probFD);
         return probFD;
         break;
      case BETA:
         if(DEBUG) printf("Beta distribution\n");
         /* -- beta distribution */
         switch(subtype)
            {
            case TAKAO:
               a1 = 0.7; a2 = -0.87; a3 = 0.0; a4 = 0.0; b1 = 2.30; b2 = -3.84; b3 = 0.0; b4 = 0.0;
               break;
            case JNES:
               a1 = 0.7; a2 = -0.81; a3 = -1.25; a4 = 0.0; b1 = 2.10; b2 = -3.84; b3 = 1.0; b4 = 0.0;
               break;
            case MOSS:
               a1 = 0.713; a2 = 0.901; a3 = 0.0; a4 = 0.0; b1 = 1.74; b2 = -1.86; b3 = 0.0; b4 = 0.0;
               if(DEBUG) printf("Moss beta\n");
               break;
            case YOUNGS:
               a1 = -0.705; a2 = 1.138; a3 = 0; a4 = 0; b1 = 0.421; b2 = -0.257; b3 = 0; b4 = 0;
               if(DEBUG) printf("Youngs beta\n");
               break;
            default:
               printf("FD sub relation type not known or incompatible: %d\n",subtype);
               return -1.;
            }
         a=exp(a1 + a2 * loverL + a3 * pow(loverL,2) + a4 * pow(loverL,3));
         b=exp(b1 + b2 * loverL + b3 * pow(loverL,2) + b4 * pow(loverL,3));
         *logD=log(a/(a+b))+logD_max;
         if(logTarget > logD_max) 
            {
            probFD=0.0;
            }
         else
            {
            x=exp(logTarget-logD_max);
            probFD = 1.-betai(a,b,x);
            }
         *Disp=exp(*logD);
         return probFD;
         break;
      case GAMMA:
         if(DEBUG) printf("Gamma distribution\n");
         /* -- gamma distribution */
         switch(subtype)
            {
            case TAKAO:
               a1 = 0.7; a2 = 0.34; a3 = a4 = 0.0; b1 = -1.40; b2 =  1.82; b3 = b4 = 0.0;
               break;
            case JNES:
               a1 = 0.5; a2 = 2.23; a3 = -4.71; a4 = 0.0; b1 = -1.15; b2 =  1.6; b3 = -0.15; b4 = 0.0;
               break;
            case MOSS:
               a1 = 0.574; a2 = -2.29; a3 = 19.9; a4 = -30.4; b1 = -1.05; b2 = 6.60; b3 = -34.6; b4 = 50.3;
               if(DEBUG) printf("Moss gamma\n");
               break;
            case YOUNGS:
               a1 = -0.193; a2 = 1.628; a3 = a4 = 0.0; b1 = 0.009; b2 = -0.476; b3 = b4 = 0.0;
               if(DEBUG) printf("Youngs gamma\n");
               break;
            default:
               printf("FD sub relation type not known or incompatible: %d\n",subtype);
               return -1.;
            }
         a=exp(a1 + a2 * loverL + a3 * pow(loverL,2) + a4 * pow(loverL,3));
         b=exp(b1 + b2 * loverL + b3 * pow(loverL,2) + b4 * pow(loverL,3));
         x=exp(logTarget-logD_av)/b;
         *logD=log(a*b)+logD_av;
         probFD=1.-gammp(a,x);
         break;
      default: 
         printf("FD relation type not know: %d\n",type);
         return -1.;
      }
   *Disp=exp(*logD);
   return probFD;
   }

double probSurfRup(float mag, int PSR)
   /* function to compute the probability of surface rupture */
   /* PSR refers to the functional form */
   /* -- compute probability of rupture --
      --  0 - normal (Youngs et al., 2003 - extensional)
      --  1 - mixed  (Wells and Coppersmith, 1993 - global)
      --  2 - thrust (Moss and Ross, 2011)
      --  3 - mixed  (Takao et al., 2013)
      --  4 - thrust stiff-soil (Moss et al., 2013)  
      --  5 - thrust soft-soil (Moss et al., 2013)  
      --  6 - strike-slip stiff-soil (Moss et al., 2013)
      --  7 - strike-slip soft-soil (Moss et al., 2013)
      --  9 - all type (Pizza et al., 2023)
      --  9 - normal (Pizza et al., 2023)
      -- 10 - reverse (Pizza et al., 2023)
      -- 11 - strike-slip (Pizza et al., 2023)

      Missing:
      --  * - thrust (Aliso Canyon, from depth distribution))
   */
   {
   float a[] = {-12.5300, -12.5100,  -7.3000, -32.0300, -13.9745,  -6.2548, -11.4071, -12.2908, -14.470, -13.500, -10.750, -28.560 };
   float b[] = {  1.9210,   2.0530,   1.0300,   4.9010,   2.1395,   0.8308,   1.8465,   1.9520,   2.177,   2.159,   1.427,   4.436 };
   double probSurfRup;
   probSurfRup=exp(a[PSR]+b[PSR]*mag)/(1.+exp(a[PSR]+b[PSR]*mag));
   return probSurfRup;
   }

double probFDDistr(float mag, double R,  double logD_av, double logD_mx,double logTarget, int FDM, double *Disp, int h)
   /* probability of distributed displacement exceedance from Petersen et al. (2011) */
   {
   double probFD, logR;
   float a,a1,a2,a3,a4,b,b1,b2,b3,b4,c1,c2,c3;
   enum EquType type;
   enum EquSubType subtype;
   float sigTrunc=12.;
   double x,sigD, logD,logD1;
   struct FDEq *D;
 
   D=InitEq();

   logR=log(R);
   type=D[FDM].type;
   subtype=D[FDM].subtype;
   switch(type)
      {
      case NORM:
         if(DEBUG) printf("Normal distribution\n");
         switch(subtype)
            {
            case QUAD:
            case ELLI:
            case BIL:
               if(DEBUG) printf("Bilinear\n");
               if(logR< 3.0) logR=3.0;
               logD=1.4016*mag-0.1671*logR-6.7991-log(100.);
               sigD=1.1193;
               break;
            case QUADN:
            case ELLIN:
            case BILN:
               if(DEBUG) printf("Bilinear normalized\n");
               if(logR< 3.0) logR=3.0;
               logD=logD_av-0.1826*logR-1.5471;
               sigD=1.0013;
               break;
            default:
               printf("FD distributed sub relation type not know: %d\n",subtype);
               return -1.;
            }
         probFD = pXceed ( logD, logTarget,sigD, sigTrunc);
         break;
      case PEMG:
         switch(subtype)
            {  
            case CHIOU:
               if(DEBUG) printf("Chiou et al. 2023\n");
               /*logD=logD_av-0.1826*logR-1.5471;
               sigD=1.0013;*/
               if(logR< 3.0) logR=3.0;
               logD=1.4016*mag-0.1671*logR-6.7991-log(100.);
               sigD=1.1193;
               break;
            default: 
               printf("FD sub relation type not know: %d\n",subtype);
               return -1.;
            }
         probFD = pXceed ( logD, logTarget,sigD, sigTrunc);
         break;
      case BETA:
      case GAMMA:
         if(DEBUG) printf("Youngs or Takeo distributed dispalcement\n");
         switch(subtype)
            {
            case YOUNGS:
               if(DEBUG) printf("Youngs gamma\n");
               /* for 95th percentile */
               if(h == 1) 
                  {
                  a=2.5; b1=0.35; b2=-0.091; b3=5.535;
                  logD=logD_av+0.35-0.091*R;
                  }
               else
                  {
                  a=2.5; b1=0.16; b2=-0.137; b3=5.535;
                  logD=logD_av+0.16-0.137*R;
                  }
               /* these need to be replaced */
               b=b1*exp(b2*R/1000.)/b3;
               logD1=log(a*b)+logD_av;
               x=exp(logTarget-logD_mx)/b;
               break;
            case TAKAO:
               a=1.282; b1=0.644; b2=-0.17; b3=1.000;
               logD=log(a*b)+logD_av;
               b=b1*exp(b2*R/1000.)/b3;
               x=exp(logTarget-logD_av)/b;
               break;
            default:
               printf("FD distributed sub relation type not know: %d\n",subtype);
               return -1.;
            }
         probFD=1.-gammp(a,x);
         break;
      default: 
         printf("FD distributed relation type not know: %d\n",type);
         return -1.;
      }
   *Disp=exp(logD);
   return probFD;
   }

double probSurfRupR(double R, double L, int iAcc, int iCmplx, int FDM)

   /* probability of surface rupture between distance Rmin and Rmax from the mapped fault trace */
   /* this factor accounts for mapping error                                                    */
   /* original form is the normal distribution                                                  */
   /* Updates:                                                                                  */
   /*          20240216 added double exponential following Chen                                 */
   /*                                                                                           */
   /* iCmplx: 0 - all, 1 - simple,   2 - complex                                                */
   /* iAcc:   0 - all, 1 - accurate, 2 - approximate, 3 - concealed, 4 - inferred               */

   {
   double sigma[3][5] = {{ 52.92, 26.89, 43.82, 65.52, 72.69},
                         { 52.92, 26.89, 43.82, 61.92, 49.57},
                         { 52.92, 26.89, 43.82, 116.2,116.35}};
   double sig, trunc=10.0,probR;
   double f;
   /* set GMPE types */
   struct FDEq *D;
   enum EquType type;
   enum EquSubType subtype;
   enum EquMapErr maperr;
   D=InitEq();

   maperr=D[FDM].maperr;
   subtype=D[FDM].subtype;
   switch(maperr)
      {
      case GAUSS:
         sig=sigma[iCmplx][iAcc];
         if(L>0.0)
            {
            probR=pXceed(0.0, R-L/2. ,sig, trunc)-pXceed(0.0,R+L/2., sig, trunc);
            }
         else
            {
            probR=pXceed(0.0,-R, sig, trunc)-pXceed(0.0,R, sig, trunc);;
            /*probR=normPDF(0.0,R,sig);*/
            }
         break;
      case DBLEXP:
         sig=sigma[iCmplx][iAcc];
         probR=pXceedLap(0.0, R-L/2. ,sig, trunc)-pXceed(0.0,R+L/2., sig, trunc);
         break;
      case NONE:
         probR=1.0;
         break;
      default: 
         probR=1.0;
         break;
      }
   return probR;
   }
   

double probSurfRupArea(float R,float A, double mag, int FDM, int h)
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
   double areaT[] = {2500, 10000, 62500, 250000 };
   double c1[] = { -6.988, -6.135, -4.903, -3.859 };
   double c2[] = { -1,410, -1.427, -1.459, -1.499 };
   double c3[] = { 0.2, 0.2, 0.2, 0.2 };
   double f;
   double Ml[] = {5.9, 6.2, 6.5, 6.8};
   double Mu[] = {6.6, 6.9, 7.2, 7.5};
   double Ql[] = { .50, .35, .20, .5 };
   double Qu[] = { 1., 1., 1., 1. };
   /* set GMPE types */
   struct FDEq *D;
   enum EquType type;
   enum EquSubType subtype;
   D=InitEq();

   type=D[FDM].type;
   subtype=D[FDM].subtype;
   switch(type)
      {
      case PEMG:
      case NORM:
         i=0;
         while (i<5 && area[i]< A) i++;
         if(R<r[i][0])
            {
            probSurfRup=p[i][0]/100.;
            }
         else if(R<r[i][2])
            {
            probSurfRup=interplilo(r[i],p[i],R,3)/100.;
            }
         else
            {
            probSurfRup=exp(a[i]*log(R)+b[i]);
            }
         break;
      case BETA:
      case GAMMA:
         switch(subtype)
            {
            case YOUNGS:
               f=2.06+(-4.63+0.118*mag+0.682*h)*log(R/1000.+3.32);
               break;
            case TAKAO:
               if(R > 0.0)
                  {
                  /* This is Takao's P2d term */
                  i=0;
                  while (i<4 && areaT[i]< A) i++;
                  /* f=-3.839+(-3.886+0.350*mag)*log(R/1000.+0.200); */
                  f=c1[i]+c2[i]*log(R/1000.+c3[i]);
                  }
               else
                  {
                  /* This is Takao's P2p term */
                  probSurfRup=0.0;
                  for(i=0;i<4;i++)
                     {
                     if(mag< Ml[i])
                        {
                        probSurfRup += .25 * Ql[i];
                        }
                     else if(mag< Mu[i])
                        {
                        probSurfRup += .25 * (Ql[i]+(Qu[i]-Ql[i])*(mag-Ml[i])/(Mu[i]-Ml[i]));
                        }
                     else
                        {
                        probSurfRup += .25 * Qu[i];
                        }
                     }
                  return probSurfRup;
                  }
               break;
            default:
               printf("FD distributed sub relation type not know: %d\n",subtype);
               return -1.;
            }
         probSurfRup = exp(f)/(1.+exp(f));
         break;
      default: 
         probSurfRup=1.0;
         break;
      }
   return probSurfRup;
   }
