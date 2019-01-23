/*---------------------------------------------------------------

  Program for checking velocity in Locating Source of Noise 

  Author : ChenXiaohan
  
  Revision History:
    2013.04.08 : initial coding

-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"
#include "sacio.h"
#include "z_sac.h"

double distaz(double,double,double,double);

int proc( char *f1, char *f2, int nn ,float coef, float SNR,float * tt_c, float * tt_s){
  int i,j,m,n1,n2,n,nlen,nwin,wlen,nfft,p;
  float *ar1,*ar2,*tmp,*ar,*ev;
  float max_c,max_s,fac1=0,fac2=0;
  char error[256];
  SACHEAD hd1,hd2;

  ar1 = (float *)malloc(nn*sizeof(float));
  ar2 = (float *)malloc(nn*sizeof(float));
  ar1 = read_sac(f1,&hd1);
  ar2 = read_sac(f2,&hd2);
  n1 = nn;
  n2 = nn;
  n = n1 + n2 - 1;
  nlen = n1;
  nwin = 1;
  wlen = nlen;
  nfft = 0;
  tmp = (float *)malloc(4*nlen*sizeof(float));
  crscor(ar1, ar2, nlen,
         nwin, wlen, SAC_RECTANGLE,
         tmp, &nfft, error, 256);
  ar = (float *)malloc(n*sizeof(float));
  for(i=0; i<=n1-2; i++)
    ar[i] = tmp[nfft-n1+i+1];
  for(i=0; i<=n2-1; i++)
    ar[n1+i-1] = tmp[i];

  ev = (float *)malloc(n*sizeof(float));
  envelope(n,ar,ev);

  for(i=0; i<n1; i++) fac1 += ar1[i]*ar1[i];
  for(i=0; i<n2; i++) fac2 += ar2[i]*ar2[i];
  for(i=0; i<n; i++) ar[i]/=sqrt(fac1*fac2);

  max_c = -1e32;
  max_s = -1e32;
  for(i=0; i<n; i++){
    if (ar[i]>max_c){
      max_c = ar[i];
      *tt_c = (i+1-n1)*hd1.delta;
      }
    if (ev[i]>max_s){
      max_s = ev[i];
      *tt_s = (i+1-n1)*hd1.delta;
      m = i;
      }
    }
  
  if ( max_c<coef || max_s<1)
    p = 0;
  else{
    int k,s[2],t[2],count=0;
    float sum=0,d1;
    d1 = hd1.delta;
    k = (int)(m - 510/d1); s[0] = (k<0)?0:k;
    k = (int)(m - 10 /d1); t[0] = (k<0)?0:k;
    k = (int)(m + 10 /d1); s[1] = (k>=n)?(n-1):k;
    k = (int)(m + 510/d1); t[1] = (k>=n)?(n-1):k;
    for(i=0; i<2; i++)
    for(j=s[i]; j<t[i]; j++){
      sum += ev[j];
      count += 1;
      }
    if (sum/count < SNR)
      p = 0;
    else 
      p = 1;
    } 
  free(tmp);
  free(ar);
  free(ar1);
  free(ar2);
  free(ev);
  return(p);
}

int main(int argc, char **argv){
  int i,k,m,NUM,n,NUM_block;
  int err,D;
  float x1,y1,x2,y2,stlo[400],stla[400],dist[400];
  float coef, SNR, tt_c, tt_s;
  char fname[400][128],f1[128],f2[128],stnm[400][128];
  FILE *fin,*fout,*fout2,*fout3,*fout4;
  SACHEAD hd;

  fin = fopen("result_rdseed","r");
  fscanf(fin,"%d",&NUM);
  for(i=0; i<NUM; i++){
    fscanf(fin,"%s",fname[i]);
    read_sachead(fname[i],&hd);
    stlo[i] = hd.stlo;
    stla[i] = hd.stla;
    dist[i] = hd.dist;
    n = hd.npts;
    strncpy(stnm[i],hd.kstnm,4);
    stnm[i][4] = '\0';
    }
  fclose(fin);

  for(i=1; i<argc; i++)
    if (argv[i][0] == '-')
      switch(argv[i][1]){
	case 'R':
	  sscanf(&argv[i][2],"%s",f2);
	break;
	case 'C':
	  sscanf(&argv[i][2],"%f",&coef);
	break;
	case 'S':
	  sscanf(&argv[i][2],"%f",&SNR);
	break;
	case 'D':
	  sscanf(&argv[i][2],"%d",&D);
	break;
	}

  fout4 = fopen("report_check_veloity","w");
  fin = fopen(f2,"r");
  fscanf(fin,"%d",&NUM_block);
  for(i=1;i<=NUM_block;i++){
    fscanf(fin,"%f %f %f %f",&x1,&x2,&y1,&y2);
    sprintf(f1,"tt_dist_%d",i);
    fout = fopen(f1,"w");
    sprintf(f1,"label_%d",i);
    fout2 = fopen(f1,"w");
    sprintf(f1,"pairs_%d",i);
    fout3 = fopen(f1,"w");
    for(k=0;k<NUM-1;k++)
    for(m=k+1;m<NUM;m++)
      if ( stlo[k]>=x1 && stlo[k]<=x2 )
      if ( stla[k]>=y1 && stla[k]<=y2 )
      if ( stlo[m]>=x1 && stlo[m]<=x2 )
      if ( stla[m]>=y1 && stla[m]<=y2 ){
	err = proc(fname[k],fname[m],n,coef,SNR,&tt_c,&tt_s);
	if ( err == 0 ) continue;

	float dist_display;
	if ( D==0 )
	  dist_display = fabs(dist[k]-dist[m]);
	else
	  dist_display = (float)distaz( (double)stla[k],
					(double)stlo[k],
					(double)stla[m],
					(double)stlo[m] );

	fprintf( fout,"%f %f\n",fabs(tt_c),dist_display );
	sprintf( f1,"%s-%s",stnm[k],stnm[m]);
	fprintf( fout2,"%f %f 3 0 4 BL %s\n",fabs(tt_c),dist_display,f1);
	fprintf( fout3,">\n%f %f\n%f %f\n",stlo[k],stla[k],stlo[m],stla[m]);
	fprintf( fout4,"%s %.2f\n",f1,tt_c-tt_s);

	}
    fclose(fout);
    fclose(fout2);
    fclose(fout3);
    }
  fclose(fin);  
  fclose(fout4);

  return 0;
}
