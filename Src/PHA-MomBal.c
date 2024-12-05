      subroutine Set_Rates ( nParamVar, magRecur, rate, beta, minMag,
     1           maxMag, iFlt, faultWidth, nWidth, fLength, slipRate,
     2           recP1, recP2, recP3, recP4,recP5,magStep )
      include 'pfrisk.h'
c     implicit none
      real magRecur(MAX_FLT,MAXPARAM,MAX_WIDTH),
     1     beta(MAX_FLT,MAXPARAM,MAX_WIDTH),
     1     minMag(MAX_FLT), maxMag(MAX_FLT,MAXPARAM,MAX_WIDTH),
     1     rate(MAXPARAM,MAX_WIDTH),
     2     faultWidth(MAX_FLT,MAX_WIDTH), fLength,
     2     slipRate(MAX_FLT,MAXPARAM,MAX_WIDTH),
     3     recP1(MAX_FLT,MAXPARAM,MAX_WIDTH),
     3     recP2(MAX_FLT,MAXPARAM,MAX_WIDTH),
     4     recP3(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real recP4(MAX_FLT,MAXPARAM,MAX_WIDTH),
     4     recP5(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real magStep(MAX_FLT)
      real mU, mL, momentMU, rigidity, momentRate1, momentRate2
      real b, c, beta1, t1, t2, t3, t4, t5, c1, deltaM1, deltaM2
      real mean, sigma, zmagL, zmagU, pmagL, pmagU, dd, mag
      real mU1, mL1, sum
      integer iParam, iFlt, i, nParamVar(MAX_FLT,1), nWidth(1)
      integer nmstep

double MomBalMX(FltArea,SlipRate,minMag,maxMag,cMag,dMag,Sigma)
   /* calculations based on slip rate constraint */
   /* Max mag model                              */
   {
   double Mu = 3.0e11;
   double rate;

   MomRate=FltArea*SlipRate*Mu*1.e9;

   nstep = (maxMag - minMag)/dMag;
   sum = 0.;
   Mag = minMag +dMag/2.;
   for(i = 0;i<nstep;i++)
      {
      sum+=probMagMX(Mag,cMag,minMag,maxMag,dMag, Sigma)* pow(10.,1.5*Mag+16.1);
      Mag += dMag;
      }
   rate = MomRate/sum;
   return rate;
   }

double MomBalGR(FltArea,SlipRate,minMag,maxMag,cMag,dMag,Beta)
   {
   double momentRate, momentRate1,c1,t1,t2,t3;
   MomRate=FltArea*SlipRate*Mu*1.e9;
   c1 = 1.5*alog(10.);
   t1 = beta * exp(beta*minMag+c1*10.7) /(1. - exp(-beta*(cMag-minMag)));
   t2 = -Beta+c1;
   t3 = ( exp(t2*maxMag) - exp(t2*minMag) );
   momentRate1 = t1 / t2 * t3;
   rate = momentRate / momentRate1;
   }

double MomBalCH(double FltArea,double SlipRate,double minMagCH, double minMag,double maxMag,double cMag,double dMag,double Beta)
   {
   double momentRate, momentRate1,c1,t1,t2,t3;
   deltaM1 = recP1(iFlt,iParam,i)
   deltaM2 = recP2(iFlt,iParam,i)
   beta1 = beta(iFlt,iParam,i)
   c1 = 1.5*alog(10.);
   t1 = exp(-beta1*(maxMag-deltaM1-minMag))
   t2 = exp(c1*(mL+10.7))
   t3 = exp(c1*(mU+10.7))
   t4 = exp(c1*(mU-deltaM1+10.7))
   t5 = exp(-beta1*(maxMag-deltaM1-deltaM2-minMag))
   momentRate1 = (beta/(c1-beta) * (t1*t4-t2)+ beta*t5/c1*(t3-t4) )/ ( (1-t1) + beta1*t5*deltaM1 )
   rate = momentRate / momentRate1
   }
   
   
double SetRate()
   {
   if(rate < 0.0)
      {
      if(GR || MX)
         {
         rate=-1.*slipRate;
      }
      else if (CH)
         {
         if ( minMag .ge. minMag) 
            {
            rate(iParam,i) = -1./slipRate(iFlt,iParam,i)
            }
         else
            {
            x = exp(beta1*(maxMag-deltaM1-deltaM2-minMag))/deltaM1
            rate = -1./slipRate *(1 + x/beta*(1-exp(-beta*(maxMag-deltaM1-minMagL))))
         }
      }
   else
      {
      switch(MagDisTyp)
      }
   }
