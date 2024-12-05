#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include "PHA.h"

#define ITMAX 100
#define NR_END 1
#define FREE_ARG char*
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

/*int main()
   {
   int i;
   double sigma=5 , mean=0.0, L=200.;
   double logTarget,p,p1,p2,w;
   for(i=-60;i<=60;i++)
      {
      logTarget=i*2;
      w = (logTarget-mean)/sigma;
      p = pXceed(mean,logTarget,sigma,12.);
      p1= normPDF(mean,logTarget,sigma);
      p2=pXceed(0.0,logTarget-L/2., sigma, 12.)-pXceed(0.0,logTarget+L/2, sigma, 12.);
      printf("%12.9e %12.9e %12.9e %12.9e\n",logTarget,p,p1,p2);
      }
   }*/

double interp(double *x,double *y,double x0,int nx)
   {
/*    linear interpolation/extrapolation. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996 */
   int i,j,isign;
   double interp;
   j=0;
   isign=1;
   if(x[0] >  x[nx-1]) isign=-1;
   while(isign*(x0-x[j]) > 0.0 ) j++;
   j--;
   /*printf("j=%5d %12.6e %12.6e\n",j,x[j],x[j+1]);*/

   if(j >=  nx-1) 
      {
      interp=y[nx-1]+(y[nx-1]-y[nx-2])*(x0-x[nx-1])/(x[nx-1]-x[nx-2]);
      }
   else
      {
      if(j <  0) j=0;
      interp=y[j]+(y[j+1]-y[j])*(x0-x[j])/(x[j+1]-x[j]);
      }
   return interp;
   }

double interpl(double *xx,double *yy,double xx0, int nx)
   {
   /*    logarithmic interpolation of function y(x(i))
         where both x should be mononous and both x and y should
         be positive. datapoint outside the interval x(1)-x(nx)
         are extrapolated, using the nearest dy/dx.
         
                                     Hong Kie Thio, January 1996*/

   double x0,*x,*y,interpl;
   int i;

   x = (double *) calloc(nx, sizeof(double));
   y = (double *) calloc(nx, sizeof(double));

   for(i=0;i<nx;i++)
      {
      if(xx[i] <  .0  ||  yy[i] <   .0) 
         {
         printf("%f %f Negative x or y, cannot take log\n",xx[i],yy[i]);
         return -1;
         }
      else
         {
         if(yy[i] == 0.0) yy[i]=1.e-308;
         if(xx[i] == 0.0) xx[i]=1.e-308;
         x[i]=log10(xx[i]);
         y[i]=log10(yy[i]);
         }
      }
   x0=log10(xx0);
   interpl=interp(x,y,x0,nx);
   if(interpl <  -308.)
      {
      interpl=-308.;
      return 1.e99;
      /*printf("Underflow, set to 1e-308\n");*/
      }
   else if(interpl >  36)
      {
      interpl=308.;
      /*printf("Overflow, set to 1e308\n");*/
      return 1.e-99;
      }
   interpl=pow(10.,interpl);
   return interpl;
   }

double interplilo(double *x,double *yy,double x0, int nx)
   {
   /*    logarithmic interpolation of function y(x(i))
         where both x should be mononous and both x and y should
         be positive. datapoint outside the interval x(1)-x(nx)
         are extrapolated, using the nearest dy/dx.
         
                                     Hong Kie Thio, January 1996*/

   double *y,interpl;
   int i;

   y = (double *) calloc(nx, sizeof(double));

   for(i=0;i<nx;i++)
      {
      if( yy[i] <   .0) 
         {
         printf("%f Negative y, cannot take log\n",yy[i]);
         return -1;
         }
      else
         {
         if(yy[i] == 0.0) yy[i]=1.e-308;
         y[i]=log10(yy[i]);
         }
      }
   interpl=interp(x,y,x0,nx);
   if(interpl <  -308.)
      {
      interpl=-308.;
      /*printf("Underflow, set to 1e-308\n");*/
      }
   else if(interpl >  36)
      {
      interpl=308.;
      /*printf("Overflow, set to 1e308\n");*/
      }
   interpl=pow(10.,interpl);
   return interpl;
   }

double pXceedLap(double Mean, double Target, double Sigma, double SigTrunc)
   {
   double x, p, pTrunc,  pXceed;
   x = (Target-Mean)/Sigma;
   p = dblexpCDF(x);
   if (x > SigTrunc)
      {
      pXceed = 0.0;
      }
   else
      {
      pTrunc =  dblexpCDF(SigTrunc);
      pXceed = (p-pTrunc)*(1.+pTrunc);
      }
   return pXceed;
   
   }


double pXceed(double Mean, double Target, double Sigma, double SigTrunc)
   /* cumlative probability for the normal distribution */
   {
   double x, p, pTrunc,  pXceed;
   x = (Target-Mean)/Sigma;
   p = normCDF(x);
   if (x > SigTrunc)
      {
      pXceed = 0.0;
      }
   else
      {
      pTrunc =  normCDF(SigTrunc);
      pXceed = (p-pTrunc)*(1.+pTrunc);
      }
   return pXceed;
   }

/* following routines from Numerical Recipes */
float gammp(float a, float x)
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}


float betai(float a, float b, float x)
{
	float betacf(float a, float b, float x);
	float gammln(float xx);
	void nrerror(char error_text[]);
	float bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double normPDF(double Mean, double Target, double Sigma)
   {
   double pNorm;
   pNorm=exp(-0.5*pow((Target-Mean)/Sigma,2))/(Sigma*sqrt(2*M_PI));
   return pNorm;
   }

double normCDF(double x)
   {
   double pNorm;
   pNorm=(1.-erf(x/sqrt(2.0)))/2.;
   return pNorm;
   }

double dblexpCDF(double x)
   {
   double pDblExp;
   if(x > 0)
      {
      pDblExp=0.5*exp(-x);
      }
   else
      {
      pDblExp=1.-0.5*exp(x);
      }
   return pDblExp;
   }


/* Numerical Recipes functions that are not directly accessed */
void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}


void gcf(float *gammcf, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

float betacf(float a, float b, float x)
{
	void nrerror(char error_text[]);
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


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

