#include <math.h>
#include <stdbool.h>

double pEMG (double x, double mu, double sigma, double nu)
   {

   /* This function computes the CDF of Exponentially-Modified-Gaussian distribuiton. */

   double  z, p, tmp;
   double alnorm(double, bool);

   z = x - mu - (pow(sigma,2) / nu);

   if (nu >   0.05 * sigma)
      {
      tmp = ( pow(mu + pow(sigma,2)/nu,2) - pow(mu,2) - 2 * x * pow(sigma,2) / nu) /( 2 * pow(sigma,2) );
      p = alnorm( (x - mu) / sigma, 0 ) - alnorm( z / sigma, 0 ) * exp(tmp);
      }
   else
      {
      p = alnorm( (x - mu) / sigma, 0 );
      }
   return p; 
   }

double alnorm ( double x, bool upper )
   {
   /*c
   cc ALNORM computes the cumulative density of the standard normal
   c  distribution.
   c
   c  Modified:
   c
   c    28 March 1999
   c
   c  Author:
   c
   c    David Hill
   c    Modifications by John Burkardt
   c
   c  Reference:
   c
   c    David Hill,
   c    Algorithm AS 66:
   c    The Normal Integral,
   c    Applied Statistics,
   c    Volume 22, Number 3, 1973, pages 424-427.
   c
   c  Parameters:
   c
   c    Input, double precision X, is one endpoint of the semi-infinite interval
   c    over which the integration takes place.
   c
   c    Input, logical UPPER, determines whether the upper or lower
   c    interval is to be integrated:
   c    1 => integrate from X to + Infinity;
   c   -1 => integrate from - Infinity to X.
   c
   c    Output, double precision ALNORM, the integral of the standard normal
   c    distribution over the desired interval.
   c*/

   double alnorm;
   double a1    =   5.75885480458  ;
   double a2    =   2.62433121679  ;
   double a3    =   5.92885724438  ;
   double b1    = -29.8213557807   ;
   double b2    =  48.6959930692   ;
   double c1    =  -0.000000038052 ;
   double c2    =   0.000398064794 ;
   double c3    =  -0.151679116635 ;
   double c4    =   4.8385912808   ;
   double c5    =   0.742380924027 ;
   double c6    =   3.99019417011  ;
   double con   =   1.28           ;
   double d1    =   1.00000615302  ;
   double d2    =   1.98615381364  ;
   double d3    =   5.29330324926  ;
   double d4    = -15.1508972451   ;
   double d5    =  30.789933034    ;
   double ltone =   7.0            ;
   double p     =   0.398942280444 ;
   double q     =   0.39990348504  ;
   double r     =   0.398942280385 ;
   double utzero=  18.66           ;
   bool up;
   double y, z;

   up = upper;
   z = x;

   if ( z < 0.0 )
      {
      up =  !up;
      z = - z;
      }
   if ( z > ltone && ( ( !up ) || utzero < z ) ) 
      {

      if ( up ) 
         {
         alnorm = 0.0;
         }
       else
         {
         alnorm = 1.0;
         }
      return alnorm;
      }
   y = 0.5 * z * z;

   if ( z < con )
      {
      alnorm = 0.5 - z * (p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
      }
   else
      {
      alnorm = r*exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
      }

   if ( !up ) 
      {
      alnorm = 1.0 - alnorm;
      }
   return alnorm;
   }
