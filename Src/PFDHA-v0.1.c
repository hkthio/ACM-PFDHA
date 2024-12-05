#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Sphere.h"

#define MAX_CONT 10
#define MAX_WID 3000
#define MAX_LEN 7500
 
/* Version 0.1: simple stepping along the fault trace, no depth stepping */
/*              that's sufficient for using the empirical relations      */

void spline(double *,double *,int, double, double, double *);
double splint(double *,double *,double *, int, double);
void getplane(struct geographic, struct geographic, struct geographic, struct geographic, struct geographic *, float *,float *,float *,float *);
float getrake(float,float,float);


int main(int ac, char **av)
   {
   /* internal variables */
   struct geographic e1,e2,e3;
   struct cartesian  c1,c2,c3;
   struct spherical  s1,s2,s3;
   double az,baz,gcarc,rot[3][3],irot[3][3];
   double pi;

   /* contour variables */
   struct Line *NewCont,Cont[MAX_CONT], SampleContour(struct Line,double );
   int np[50],nCont,center;

   /* sampling parameters */
   struct geographic so,sc[MAX_CONT][MAX_LEN],**sd;
   int nx,nn[3000],nw[3000],nwmax=0,iflag;
   float dx=.001,convergence;
   double xx1,yy1,cdist,xx,yy,x[3000],y[3000],y2[3000];
   double yp1=1.e30;
   double ypn=1.e30;
   float yprev,xprev,W;
   float avslip=1070.,slip;
   FILE *fp, *fp1, *fp2;

   /* input parameters */
   char cont_file[80];
   float conv_az, dpdq1,daz;
   float dL,dLL,dW,dmax;

   /* extrapolation */
   float hdist, dmax3;
   /* output parameters */
   struct geographic st, su[4], sf,sfc[4];
   float strike, strike0, dip, rake, L, area;

   int i,j,k,l;
   char cdum[80];
   pi=acos(-1.0);
   sd = (struct geographic **) calloc(MAX_WID,sizeof( struct geopgraphic *));
   for(k=0;k<MAX_WID;k++)
      {
      sd[k] = (struct geographic *) calloc (MAX_LEN,sizeof(struct geographic));
      }

   /* input parameters */
   printf("Contour/trace file?\n");
   scanf("%s",cont_file);
   printf("Convergence azimuth?\n");
   scanf("%f",&conv_az);
   printf("Sampling length and width?\n");
   scanf("%f %f",&dL, &dW);
   dL=dL/111.;
   printf("Seismogenic depth?\n");
   scanf("%f",&dmax);
   printf("Initial dip?\n");
   scanf("%f",&daz);
   if(daz <= 90. && daz >= 0.) dpdq1 = tan(daz*pi/180.);

   fp=fopen(cont_file,"r");
   k=0;
   fscanf(fp,"%d%*[^\n]\n",&nCont);
   NewCont=(struct Line *) calloc(nCont,sizeof(struct Line));
   for(k=0;k<nCont;k++)
      {
      fscanf(fp,"%s %d%*[^\n]\n",cdum,&Cont[k].np); 
      printf("%d\n",Cont[k].np);
      Cont[k].pt = (struct geographic *) calloc(Cont[k].np,sizeof (struct geographic));
      printf("Cont[%d] initialized\n",k);
      for(i=0;i<Cont[k].np;i++) 
         {
         fscanf(fp,"%lf %lf %lf",&Cont[k].pt[i].lon,&Cont[k].pt[i].lat,&Cont[k].pt[i].dep);
         printf("%f %f %f\n",Cont[k].pt[i].lon,Cont[k].pt[i].lat,Cont[k].pt[i].dep);
         if(Cont[k].pt[i].dep < 0.0) Cont[k].pt[i].dep *= -1;
         }
      }
   fclose(fp);
   /* resample the contour using straight interpolation */
   for(k=0;k<nCont;k++)
      {
      NewCont[k]=SampleContour(Cont[0], .01);
      printf("> %3d\n",NewCont[k].np);
      for(i=0;i<NewCont[k].np;i++)
         {
         printf("%10.5f %9.5f %7.2f\n",NewCont[k].pt[i].lon,NewCont[k].pt[i].lat,NewCont[k].pt[i].dep);
         }
      }
   }

struct Line SampleContour(struct Line Cont, double dr)
   {
   struct Line ContNew;
   int i,j,k,n_new,*n,ntot=0;
   double az,baz,gcarc, *r,rr;

   /* find number of elements for new contour */
   printf("Resample dr=%6.3f\n",dr);
   r=(double *) calloc(Cont.np,sizeof(double));
   n=(int *) calloc(Cont.np,sizeof(int));
   for(i=0;i<Cont.np-1;i++)
      {
      distaz(Cont.pt[i],Cont.pt[i+1],(r+i),&az,&baz);
      n[i]=0.5+r[i]/dr-1;
      printf("%f %f %f %d %f\n",Cont.pt[i].lon,Cont.pt[i].lat,Cont.pt[i].dep,n[i],r[i]);
      ntot+=n[i];
      }
   n[Cont.np-1]=1;
   ntot++;
   printf("%d\n",ntot);
   ContNew.np=ntot;
   ContNew.pt = (struct geographic *) calloc(ntot,sizeof (struct geographic));
   k=0;
   /* resample the original contour */
   for(i=0;i<Cont.np-1;i++)
      {
      ContNew.pt[k]=Cont.pt[i];
      k++;
      for(j=0;j<n[i]-1;j++)
         {
         distaz(ContNew.pt[k-1],Cont.pt[i+1],&gcarc,&az,&baz);
         rr=gcarc/(n[i]-j);
         ContNew.pt[k]=newloc(ContNew.pt[k-1],az,rr);
         printf("%f %f %f %d %d %d\n",ContNew.pt[k].lon,ContNew.pt[k].lat,ContNew.pt[k].dep,i,j,k);
         k++;
         }
      }
   ContNew.pt[k]=Cont.pt[Cont.np-1];
   return ContNew;
   }

#if 0
   fp=fopen("o_extra","w");
   dmax3=dmax/3.;
   hdist=dmax3/tan(daz*pi/180.);
   /* project downdip to create more contours */
   printf("hdist: %7.2f %7.2f %d\n",hdist,dmax3,ncont);
   fprintf(fp,"> %d\n",np[0]);
   for(i=0;i<np[0];i++)
      {
      fprintf(fp,"%f %f %f\n",co[0][i].lon,co[0][i].lat,co[0][i].dep);
      }
   distaz(co[0][0],co[0][np[0]-1],&gcarc,&az,&baz);
   for(k=ncont;k <(ncont+5);k++)
      {
      np[k]=np[0];
      fprintf(fp,"> %d\n",np[k]);
      for(i=0;i<np[k];i++)
         {
         /*if(i > 0 && i < np[k]-1)
            {
            distaz(co[k-1][i-1],co[k-1][i+1],&gcarc,&az,&baz);
            }
         else if (i == 0)
            {
            distaz(co[k-1][i],co[k-1][i+2],&gcarc,&az,&baz);
            }
         else
            {
            distaz(co[k-1][i-2],co[k-1][i],&gcarc,&az,&baz);
            }*/
         co[k][i]=newloc(co[k-1][i],az+90.,hdist/111.);
         co[k][i].dep=k*dmax3;
         fprintf(fp,"%f %f %f %f %d %d\n",co[k][i].lon,co[k][i].lat,co[k][i].dep,az,k,i);
         }
      }
   ncont=5;
   fclose(fp);
   return np_new;
   }

   center=ncont/2;
   center=0;
   printf("# %d contours read, center of projection is contour %d\n",ncont,center);

   /* set new coordinate system */
   e1.lat=co[center][0].lat;
   e1.lon=co[center][0].lon;
   e1.dep=0.0;
   printf("# %g %g %g\n",e1.lat,e1.lon,az);
   
   e2.lat=co[center][np[center]-1].lat;
   e2.lon=co[center][np[center]-1].lon;
   e2.dep=0.0;
   distaz(e1,e2,&gcarc,&az,&baz);
   if(conv_az > 360 || conv_az < -360.)
      {
      convergence=az-90.;
      }
   else
      {
      convergence=conv_az;
      }
   printf("# %g %g %g %g\n",e2.lat,e2.lon,az,convergence);

   /* set new basis to e1, and determine rotation matrix */
   c1=sphtocart(geotosph(e1));
   c2=sphtocart(geotosph(newloc(e1,az,90.)));
   c3=sphtocart(geotosph(newloc(e1,az+90.,90.)));
   rotmat(c1, c2, c3, rot, irot);

   fp=fopen("o_cont_r.xy","w");
   for(k=0;k<ncont;k++)
      {
      fprintf(fp,"> %3d\n",np[k]);
      for(i=0;i<np[k];i++)
         {
         cr[k][i]=rotgeo(co[k][i],rot);
         fprintf(fp,"%9.4f %8.4f  ",cr[k][i].lon,cr[k][i].lat);
         fprintf(fp,"%9.4f %8.4f %5.2f\n",co[k][i].lon,co[k][i].lat,co[k][i].dep);  
         }
      }
   fclose(fp);
   /* set number of points along length */
   fp=fopen("o_cont","w");
   for(k=0;k<ncont;k++)
      {
      /* setup the spline for this contour */
      fprintf(fp,">%d\n",k);
      for(i=0;i<np[k];i++)
         {
         x[i]=cr[k][i].lon;
         y[i]=cr[k][i].lat;
         }
      spline(x,y,np[k],yp1,ypn,y2);
      /* now do very dense sampling of spline */
      /*first compute distance along spline */
      i=0;
      j=0;
      xx=cr[k][0].lon;
      yy=cr[k][0].lat;
      cdist=0;
      while(xx<cr[k][np[k]-1].lon)
         {
         xx1=cr[k][0].lon+i*dL/1000.;;
         yy1=splint(x,y,y2,np[k],xx);
         cdist+=sqrt(pow((xx1-xx),2)+pow((yy1-yy),2));
         xx=xx1;
         yy=yy1;
         i++;
         }
      if(k==0) 
         {
         nx=cdist/dL+1;
         if(nx >= MAX_LEN)
            {
            printf("Increase MAX_LEN: %d\n",nx);
            return -1;
            }
         }
      dLL=cdist/(nx-1);
      printf("Contour %d: %d points, %7,3f %12.3f\n",k,nx,dLL,cdist);
      j=1;
      xx=cr[k][0].lon;
      yy=cr[k][0].lat;
      cdist=0;
      sc[k][0].lon=xx;
      sc[k][0].lat=yy;
      sc[k][0].dep=cr[k][0].dep;
      so=rotgeo(sc[k][0],irot);
      fprintf(fp,"%9.4f %8.4f %5.1f\n",so.lon,so.lat,so.dep);
      i=1;
      while(i < nx)
         {
         xx1=cr[k][0].lon+j*dL/1000.;;
         yy1=splint(x,y,y2,np[k],xx);
         cdist+=sqrt(pow((xx1-xx),2)+pow((yy1-yy),2));
         if(cdist > dLL)
            {
            sc[k][i].lon=xx1;
            sc[k][i].lat=yy1;
            sc[k][i].dep=cr[k][0].dep;
            so=rotgeo(sc[k][i],irot);
            fprintf(fp,"%9.4f %8.4f %5.1f\n",so.lon,so.lat,so.dep);
            cdist-=dLL;
            i++;
            }
         xx=xx1;
         yy=yy1;
         j++;
         }
      /*for(i=0;i<nx;i++)
         {
         yy=splint(x,y,y2,np[k],xx);
         sc[k][i].lon=xx;
         sc[k][i].lat=yy;
         sc[k][i].dep=cr[k][0].dep;
         so=rotgeo(sc[k][i],irot);
         fprintf(fp,"%9.4f %8.4f %5.1f\n",so.lon,so.lat,so.dep);
         xx+=dLL;
         }*/
      }
   fclose(fp);
   /* now find intermediate depth points for every nx points */
   /* since we have sampled all contours with the same x values, there is no need for horizontal spline */
   for(i=0;i<nx;i++)
      {
      for(k=0;k<ncont;k++)
         {
         distaz(sc[0][i],sc[k][i],&gcarc,&az,&baz);
         x[k]=gcarc*111.;
         y[k]=sc[k][i].dep;
         }
      spline(x,y,ncont,yp1,ypn,y2);
      j=0;
      xprev=x[0];
      yprev=y[0];
      xx=0.0;
      sd[0][i]=sc[0][i];
      printf("y[0]=%f\n",yprev);
      while(yprev< dmax)
         {
         W=0.0;
         while(W < dW)
            {
            yy=splint(x,y,y2,ncont,xx);
            W=sqrt(pow(xx-xprev,2)+pow(yy-yprev,2));
            xx=xx+dx;
            }
         j++;
         /*printf("%d %f %f %f %f %f\n",j,xx,yy,xprev,yprev,dmax);*/
         if(j >= MAX_WID)
            {
            printf("Increase MAX_WID: %d %d %f %f %f %f\n",j,MAX_WID,W,dW,yprev,dmax);
            return -1;
            }
         xprev=xx;
         yprev=yy;
         sd[j][i]=newloc(sc[0][i],az,xx/111.);
         sd[j][i].dep=yy;
         /*printf("# %d %d %g %g\n",j,i,sd[j][i].lon,sd[j][i].dep);*/
         }
      nw[i]=j;
      }
   /* now create the subfaults */
   fp=fopen("i_invall","w");
   fp1=fopen("o_subfaults","w");
   for(i=0;i<nx-1;i++)
      {
      if(nw[i]>nwmax) nwmax=nw[i];
      if(nw[i] > nw[i+1])
         {
         nn[i]=nw[i+1];
         }
      else
         {
         nn[i]=nw[i];
         }
      printf("# %d %d %d\n",i,nn[i],nw[i]);
      for(k=0;k<nn[i]-1;k++)
         {
         sfc[0]=sd[k][i];
         sfc[1]=sd[k][i+1];
         sfc[2]=sd[k+1][i+1];
         sfc[3]=sd[k+1][i];
         sf.lon=(sfc[0].lon+sfc[1].lon+sfc[2].lon+sfc[3].lon)/4.;
         sf.lat=(sfc[0].lat+sfc[1].lat+sfc[2].lat+sfc[3].lat)/4.;
         sf.lat=(sfc[0].dep+sfc[1].dep+sfc[2].dep+sfc[3].dep)/4.;
         printf(">%d\n",k);
         fprintf(fp1,"> %3d %3d\n",i,k);
         for(j=0;j<4;j++) 
            {
            printf("%9.4f %8.4f\n",sfc[j].lon,sfc[j].lat);
            su[j]=rotgeo(sfc[j],irot);
            fprintf(fp1,"%9.4f %8.4f\n",su[j].lon,su[j].lat);
            }
         getplane(su[0],su[1],su[2],su[3],&st,&strike,&dip,&L, &area);
         if(strike > 360) strike -=360.;
         if(strike < -180) strike +=360.;
         if(k==0) strike0=strike;
         if(st.lon < 0.0) st.lon += 360.;
         convergence=strike0-90.;
         rake=getrake(strike, dip,convergence);
         slip=100.;
         fprintf(fp,"%9.4f %8.4f %7.2f %6.1f %5.1f %6.1f %7.2f %7.2f %7.2f %3d %3d 1\n",st.lon,st.lat,st.dep,strike,dip,rake,slip,L,area/L,i,k);
         }
      }
   fp2=fopen("o_rates","w");
   fprintf(fp2,"%4d %4d\n",nx-1,nwmax-1);
   for(i=0;i<nx-1;i++)
      {
      for(k=0;k<nwmax-1;k++)
         {
         iflag=0;
         if(k<nw[i]) iflag=1;
         fprintf(fp2,"%d ",iflag);
         }
      fprintf(fp2,"\n");
      }
   fclose(fp);
   fclose(fp1);
   fclose(fp2);
   }

/*    -- Determine best fitting plane through four corner points --
      ------------------------------------------------------------- */
void getplane(struct geographic g1, struct geographic g2, struct geographic g3, struct geographic g4, struct geographic *ctr, float *strike,float *dip,float *length,float *area)
   {
   float sgn;
   int i;
   struct cartesian s[4],cr = {0.0,0.0,0.0};
   double conv,pi;
   double r,colat,rlon;
   double nx,ny,nz,ox,oy,oz,px,py,pz,qx,qy,qz,rx,ry,rz,sx,sy,sz,sr;
   double tx,ty,tz;
   double area1, area2;
   double ax,ay,az,bx,by,bz;
   double l1,l2,l3,l4;
   float strk;

   pi=acos(-1.0);
   conv=pi/180.;

   /*  - change to cartesian coordinates - */
   s[0]=sphtocart(geotosph(g1));     
   s[1]=sphtocart(geotosph(g2));     
   s[2]=sphtocart(geotosph(g3));     
   s[3]=sphtocart(geotosph(g4));     

    /* - centroid of plane is also the normal of the earth at that location - */
    for(i=0;i<4;i++)
       {
       cr.x+=0.25*s[i].x;  
       cr.y+=0.25*s[i].y;  
       cr.z+=0.25*s[i].z;  
       }
   rx=cr.x;
   ry=cr.y;
   rz=cr.z;
   /* -- center of plane - */
   *ctr=sphtogeo(cartosph(cr));

   ax=s[1].x-s[0].x;
   ay=s[1].y-s[0].y;
   az=s[1].z-s[0].z;
   bx=s[2].x-s[0].x;
   by=s[2].y-s[0].y;
   bz=s[2].z-s[0].z;
   ox=ay*bz-az*by;
   oy=az*bx-ax*bz;
   oz=ax*by-ay*bx;
   area1=0.5*sqrt(ox*ox+oy*oy+oz*oz);
        
   ax=s[3].x-s[2].x;
   ay=s[3].y-s[2].y;
   az=s[3].z-s[2].z;
   bx=s[0].x-s[2].x;
   by=s[0].y-s[2].y;
   bz=s[0].z-s[2].z;
   px=ay*bz-az*by;
   py=az*bx-ax*bz;
   pz=ax*by-ay*bx;
   area2=0.5*sqrt(px*px+py*py+pz*pz);

   /* - cumulative values - */
   *area=area1+area2;
   qx=px+ox;
   qy=py+oy;
   qz=pz+oz;
        
   /* - get dip - */
   *dip=acos((rx*qx+ry*qy+rz*qz)/sqrt((rx*rx+ry*ry+rz*rz)*(qx*qx+qy*qy+qz*qz)));
   *dip=*dip/conv;
   if(*dip > 90.) *dip=180.-*dip;

   /* - get strike - */
   /* - first get local north vector - */
   if(rz == 0.0) 
      {
      nz=1.0;
      nx=0.0;
      ny=0.0;
      }
   else if(rz > 0.0) 
      {
      nx=-rx;
      ny=-ry;
      nz=(rx*rx+ry*ry)/rz;
      }
   else
      {
      nx=rx;
      ny=ry;
      nz=-(rx*rx+ry*ry)/rz;
      }
   /* - get the strike vector - */
   sx=(ry*qz-rz*qy);
   sy=(rz*qx-rx*qz);
   sz=(rx*qy-ry*qx);
   sr=sqrt(sx*sx+sy*sy+sz*sz);
   sx=sx/sr;
   sy=sy/sr;
   sz=sz/sr;

   /* - use inner product to get strike angle - */
   strk=acos((nx*sx+ny*sy+nz*sz)/sqrt((nx*nx+ny*ny+nz*nz)*(sx*sx+sy*sy+sz*sz)));
   /* - sign of strike - */
   tx=ny*sz-nz*sy;
   ty=nz*sx-nx*sz;
   tz=nx*sy-ny*sx;
   sgn=rx*tx+ry*ty+rz*tz;
   strk=fabs(strk)*sgn/fabs(sgn);
   *strike=-strk/conv-180.;

   /* - get length along strike - */
   l1=fabs((s[0].x-rx)*sx+(s[0].y-ry)*sy+(s[0].z-rz)*sz);
   l2=fabs((s[1].x-rx)*sx+(s[1].y-ry)*sy+(s[1].z-rz)*sz);
   l3=fabs((s[2].x-rx)*sx+(s[2].y-ry)*sy+(s[2].z-rz)*sz);
   l4=fabs((s[3].x-rx)*sx+(s[3].y-ry)*sy+(s[3].z-rz)*sz);
   *length=(l1+l2+l3+l4)/2.;
   }

/*    -------------------------------------------------------------------------- */
/*    -- Get rake angle from horizontal convergence vector and fault geometry -- */
float getrake(float str,float dip,float azi)
   {
   double sn,se,sz,dn,de,dz,nn,ne,nz,nr;
   double pn,pe,pz,qn,qe,qz,tn,te,tz;
   double dp,pi,conv;
   float sgn, rak;

   printf("#getrake: %g %g %g\n",str,dip,azi);
   /* -- local cartesian coordinates, z=up -- */
   pi=acos(-1.0);
   conv=pi/180.;
   /* convergence vector */
   pn=cos(azi*conv);
   pe=sin(azi*conv);
   pz=0.0;
   /* fault plane -- */
   /* strike direction -- */
   sn=cos(str*conv);
   se=sin(str*conv);
   sz=0.0;
   /* dip direction -- */
   dn=cos((str+90.)*conv)*cos(dip*conv);
   de=sin((str+90.)*conv)*cos(dip*conv);
   dz=-sin(dip*conv);
   /* normal to the plane (upward) -- */
   nn=de*sz-dz*se;
   ne=dz*sn-dn*sz;
   nz=dn*se-de*sn;
   /* nornalize (just in case) -- */
   nr=sqrt(nn*nn+ne*ne+nz*nz);
   nn=nn/nr;
   ne=ne/nr;
   nz=nz/nr;
   /* projection of convergence vector */
   dp=pn*nn+pe*ne+pz*nz;
   qn=pn-dp*nn;
   qe=pe-dp*ne;
   qz=pz-dp*nz;
   printf("#getrake: %g %g %g %g %g %g %g %g %g\n",pn,pe,pz,nn,ne,nz,qn,qe,qz);
   /* rake angle -- */
   rak=acos((qn*sn+qe*se+qz*sz)/sqrt((qn*qn+qe*qe+qz*qz)*(sn*sn+se*se+sz*sz)));
   /* sign of strike */
   tn=se*qz-sz*qe;
   te=sz*qn-sn*qz;
   tz=sn*qe-se*qn;
   sgn=nn*tn+ne*te+nz*tz;
   rak=fabs(rak)*sgn/fabs(sgn);
   rak/=conv;
   return rak;
   }

void reverse_d( double *vec, int nvec)
/* reverses order of vector */
   {
   int i;
   double dum;
   for(i=0;i<nvec/2;i++)
      {
      dum= vec[i];
      vec[i]=vec[nvec-1-i];
      vec[nvec-1-i]=dum;
      }
   }

#endif
