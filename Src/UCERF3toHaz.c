#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "CatUtil.h"

#define MAXLINE 2048
#define MAXSEG  3000

float getMag(float);

int main(int argc, char **argv)
   {
   struct Fault {
      int id;
      int np;
      struct FP fp;
      struct geographic *tr;
      } *flt;
   struct geographic UL, BR,x,pt;
   double dr,dist,dist1;
   int   nlon,nlat,l,iplp;

   float dip, dipdir, top, bot, rake;
   float mu;
   
   int id[3000],idx,nsub,np[3000],npts;
   char* element,coordfile[128],paramfile[128];
   float fdum[120];

   int i,j,k,*nseg,idum,nev,nevss,fmin,fmax,model;
   FILE *fp;
   char segfile[128],rupfile[128],line[MAXLINE];
   int *segid;

   flt = (struct Fault *) calloc(MAXSEG,sizeof (struct Fault ));
   UL.lon=-119.5;
   UL.lat=34.8;
   BR.lon=-117.1;
   BR.lat=33.1;
   dr=.001;
   mu=30.e9;
   model=1;
   if(model == 1) 
      {
      fmin=2301;
      fmax=2305;
      sprintf(coordfile,"FM3.1_FaultTraceCoords.csv");
      sprintf(paramfile,"FM3.1_FaultDipRakeZ.csv");
      sprintf(segfile,"FM3.1_FaultIDs.txt");
      }
   else
      {
      fmin=2362;
      fmax=2370;
      sprintf(coordfile,"FM3.2_FaultTraceCoords.csv");
      sprintf(paramfile,"FM3.2_FaultDipRakeZ.csv");
      sprintf(segfile,"FM3.2_FaultIDs.txt");
      }

   /* read segment coordinates */
   printf("Reading segment coordinates %s\n",coordfile);
   i=0;
   /* KEEP THIS */
   fp=fopen(coordfile,"r");
   while(fgets(line,MAXLINE,fp) != NULL)
      {
      element=strtok(line,",");
      sscanf(element,"%d,",&k);
      if(k>=MAXSEG)
         {
         printf("Segment index exceeds MAXSEG: %d %d\n",k,MAXSEG);
         return -1;
         }
      flt[k].id=k;
      j=0;
      while((element=strtok(NULL,",")) != NULL)
         {
         if(sscanf(element,"%f,",(fdum+j))==1) j++;
         }
      flt[k].np=j/2;
      printf("> %d %d %d\n",i,k,flt[k].np);
      flt[k].tr = (struct geographic *) calloc(flt[k].np,sizeof (struct geographic)); 
      for(j=0;j<flt[k].np;j++)
         {
         flt[k].tr[j].lon=fdum[j];
         flt[k].tr[j].lat=fdum[flt[k].np+j];
         printf("%10.5f %9.5f\n",flt[k].tr[j].lon,flt[k].tr[j].lat);
         }
      i++;
      }
   nsub=i;
   fclose(fp);
   printf("Finished reading segment coordinates: %d\n",nsub);

   fp=fopen(paramfile,"r");
   j=0;
   while(fscanf(fp,"%d, %f, %f, %f, %f, %f",&i,&dip,&dipdir,&rake,&top,&bot)==6)
      {
      if(i>MAXSEG)
         {
         printf("Segment index exceeds MAXSEG: %d %d\n",i,MAXSEG);
         return -1;
         }
      flt[i].fp.dip=dip;  
      flt[i].fp.str=dipdir-90.;  
      if(flt[i].fp.str< -180) flt[i].fp.str +=360.;
      flt[i].fp.rak=rake;
      for(k=0;k<flt[i].np;k++) flt[i].tr[k].dep=top;
      j++;
      }
   if(j != nsub)
      {
      printf("Segment and parameter files have inconsistent number of entries: %d %d\n",i,nsub);
      return -1;
      }
   fclose(fp);
   printf("Finished reading segment parameters: %d\n",nsub);

   fp=fopen("o_faults","w");
   for(i=0;i<nsub;i++)
      {
      fprintf(fp,"> %3d %3d\n",i,flt[i].np);
      for(j=0;j<flt[i].np;j++)
         {
         fprintf(fp,"%10.5f %9.5f\n",flt[i].tr[j].lon,flt[i].tr[j].lat);
         }
      }
   fclose(fp);
   
   /*
   nlon=(BR.lon-UL.lon)/dr+1;
   nlat=(UL.lat-BR.lat)/dr+1;
   fp=fopen("o_dist","w");
   for(k=0;k<nlat;k++)
      {
      pt.lat=UL.lat-k*dr;
      for(l=0;l<nlon;l++)
         {
         pt.lon=UL.lon+l*dr;
         dist=999999.0;
         for(i=0;i<nsub;i++)
            {
            for(j=0;j<flt[i].np-1;j++)
               {
               dist1 = distoline(pt,flt[i].tr[j], flt[i].tr[j+1],&x,&iplp);
               if(dist1 < dist) dist=dist1;
               }
            }
         fprintf(fp,"%10.5f %9.5f %8.3f\n",pt.lon,pt.lat,dist);
         }
      }
   fclose(fp);
   */
   }
