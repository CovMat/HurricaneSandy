/*---------------------------------------------------------------

  Program for checking cross-correlation locating Source of Noise 
 
  Author : ChenXiaohan
   
  Revision History:
    2013.04.11 : initial coding
    2013.04.13 : remove Projection code , add Cross-Correlation output
    2013.04.15 : adjust report output
    2013.04.16 : check approximate location of hurricane
    2013.04.17 : add t2 check, dist > 20km check
    2013.05.14 : (1) remove -R -V options
		             (2) hint code of (Hurricane) time stamp
		             (3) add a Hurricane Position option
		             (4) add a percentage screen output
    2013.11.04 : clean the code, for calculating the katrina data
    2013.11.05 : change the omp to inner j for loop
    2014.03.11 : (1) display the finished calculation percentage
                 (2) add filenames output
                 (3) change the lowest peak value limit to 1e-6
    2014.03.12 : (1) delete the percentage output
                 (2) output the number of pairs to a individual file
    2014.06.26 : (1) add a segment number option
                 (2) add a check option : output all cross-correlation results
    2014.07.05 : add ref. peak, coef, SNR into arguments
    2014.07.29 : little change

-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "sac.h"
#include "sacio.h"
#include "z_sac.h"
#include "cal_omp.h"

void trunc_ar(float * ar, int n){
  int i;
  for(i=0;i<n;i++){
    if ( ar[i]>100 ) ar[i]=100;
    if ( ar[i]<-100 ) ar[i]=-100;
    }
}

int main( int argc, char **argv){
  int i, NUM, count = 0, seg, check_cc = 0;
  float ref_peak = 1e-6, ref_coef = 0.2, ref_SNR = 4.5;
  FILE *fin, *fout, *fout_p;
  char fnm[500][256],label[256],f1[256];

  /* read the arguments */
  for(i=1;i<argc;i++)
    if (argv[i][0]=='-')
      switch (argv[i][1]){
        case 'L':
          sscanf(&argv[i][2],"%s",label);
        break;
        case 'T':
          sscanf(&argv[i][2],"%d",&seg);
        break;
        case 'C':
          sscanf(&argv[i][2],"%d",&check_cc);  // output cross-correlation results option, 1 for output, 0 for not [default]
        break;
        case 'R':
          sscanf(&argv[i][2],"%f/%f/%f", &ref_peak, &ref_coef, &ref_SNR); // 1e-6/0.2/4.5 [default]
        break;
        }

  /* read file names  */
  fin = fopen("result_rdseed","r");
  fscanf(fin,"%d",&NUM);
  for(i=0; i<NUM; i++){
    fscanf(fin,"%s",fnm[i]);
    sprintf(fnm[i],"%s_%d",fnm[i],seg);
    }
  fclose(fin);

  sprintf(f1,"answer_%s",label);
  fout = fopen(f1,"w");
  sprintf(f1,"pairs_%s",label);
  fout_p=fopen(f1,"w");

  /* let's do it! */
  for(i=0; i<NUM-1; i++){
    int j, n1;
    float *ar1;
    SACHEAD hd1 = sac_null;
    read_sachead(fnm[i],&hd1);
    if (hd1.t2 > -12345.) continue;
    ar1 = read_sac(fnm[i],&hd1); n1 = hd1.npts;

    trunc_ar(ar1,n1);
    
    #pragma omp parallel for firstprivate(NUM,fnm,i,n1,ar1,hd1,label,check_cc) shared(count,fout,fout_p)
    for(j=i+1; j<NUM; j++){
      int n2, n, err;
      float *ar2, *ar3, *ar4, *ar_c, *ar_ev;
      float peak, coef, SNR, tt_c, tt_s, dist;
      SACHEAD hd2 = sac_null;
      read_sachead(fnm[j],&hd2);
      if (hd2.t2 > -12345.) continue;
      ar2 = read_sac(fnm[j],&hd2); n2 = hd2.npts;
      trunc_ar(ar2,n2);

      ar3 = crscor_0(ar1, ar2, n1, n2);
      n = n1 + n2 -1;
      ar4 = norm_0(ar1, ar2, ar3, n1, n2, n);
      ar3 = enve_0(ar3,n);
      ar_c = ar4;
      ar_ev = ar3;
 
      /* check cross-correlation */
	    err = check(ar_c,ar_ev,hd1,hd2,
                  ref_peak,ref_coef,ref_SNR,
		              &peak,&coef,&SNR,
		              &tt_c, &tt_s);
	    dist = (float)distaz( (double) hd1.stla,
			                      (double) hd1.stlo,
			                      (double) hd2.stla,
			                      (double) hd2.stlo );

      /* output cross-correlation */
      if (dist > 20)
      if (check_cc == 1){
        char cc_filename[256];
        sprintf(cc_filename,"cc_%s_%d_%d",label,i,j);
        SACHEAD cc_hd = sachdr( hd1.delta, n, hd1.delta*(1-n1) );
        strncpy( cc_hd.kstnm, hd1.kstnm, 8 );
        strncpy( cc_hd.kcmpnm, hd2.kstnm, 8 );
        cc_hd.stlo = hd1.stlo;
        cc_hd.stla = hd1.stla;
        cc_hd.evlo = hd2.stlo;
        cc_hd.evla = hd2.stla;
        cc_hd.o = 0;
        cc_hd.iztype = IUNKN;
        write_sac( cc_filename, cc_hd, ar_c );
        }

      if (err == 0)
	    if (dist > 20){
        #pragma omp atomic
          count += 1;
        #pragma omp critical (file_lock)
          {
          fprintf(fout,"%f %f %f %f %f %f %d %d %f %f %f\n",tt_c,dist,hd1.stlo,hd1.stla,hd2.stlo,hd2.stla,i,j,peak,coef,SNR);
//          fprintf(fout,"%s\n%s\n",fnm[i],fnm[j]);
          fprintf(fout_p,">\n%f %f\n%f %f\n",hd1.stlo,hd1.stla,hd2.stlo,hd2.stla);
          }
      	}
      free(ar2);
      free(ar3);
      free(ar4);
      }
    free(ar1);  
    } 

//  fprintf(fout,"0 0 0 0 0 0 0 0 0 0 0\n%d",count);
  fclose(fout);
  fclose(fout_p);

  sprintf( f1, "num_picked_pairs_%s", label );
  fout = fopen( f1 ,"w");
  fprintf(fout,"%d",count);
  fclose(fout);

  return 0;
}

int check( const float *ar_c, const float *ar_ev,
           const SACHEAD hd1, const SACHEAD hd2,
	   const float l_peak, const float l_coef, const float l_SNR,
	   float *p_peak, float *p_coef, float *p_SNR, 
	   float *p_tt_c, float *p_tt_s ){
  int i,j,m,s[2],t[2],k;
  int n,n1,n2,tmax,t1,t2,count;
  float sum,noi,d;
  float max_c, max_s, tt_c, tt_s;
   
    max_c = -1e32;
    max_s = -1e32;
    sum = 0;
    count = 0;
    n1 = hd1.npts;
    n2 = hd2.npts;
    n = n1 + n2 - 1;
    d = hd1.delta;
   
    tmax = (int)( distaz( (double) hd1.stla,
                          (double) hd1.stlo,
                          (double) hd2.stla,
                          (double) hd2.stlo )/1.0/d );
    t1 = n1-1 - tmax;
    t1 = (t1<0)?0:t1;
    t2 = n1-1 + tmax;
    t2 = (t2>n)?n:t2;
    m = t1;
    for(i=t1; i<=t2; i++){
      if (ar_c[i] > max_c){ 
	      max_c = ar_c[i]; 
	      tt_c=(i+1-n1)*d; 
	      m = i;
	      }
      if (ar_ev[i] > max_s){
	      max_s = ar_ev[i];
	      tt_s=(i+1-n1)*d;
	      }
      }
  
    if ( max_c<l_coef || max_s<l_peak )
      return 1;
    else{
      k = (int)(m - 510/d); s[0] = (k<0)?0:k;
      k = (int)(m - 10 /d); t[0] = (k<0)?0:k;
      k = (int)(m + 10 /d); s[1] = (k>=n)?(n-1):k;
      k = (int)(m + 510/d); t[1] = (k>=n)?(n-1):k;
      for(i=0; i<2; i++)
      for(j=s[i]; j<t[i]; j++){
	      sum += ar_ev[j];
	      count += 1;
	      }
      noi = sum / count;
      if ( max_s/noi<l_SNR )
	      return 1;
      }

    *p_peak = max_s; *p_coef = max_c; *p_SNR = max_s/noi;
    *p_tt_c = tt_c; *p_tt_s = tt_s;
    return 0;
}


float* enve_0( float *ar,
               const int n ){
  float *ev, *swap; 
  ev = (float *)malloc(n*sizeof(float));
  envelope(n, ar, ev);
  swap = ar; ar = ev; ev = swap;
  free(ev);
  return ar;
}

float *norm_0( const float *x,
	       const float *y,
	       const float *z,
	       const int   n1,
	       const int   n2,
	       const int   n){
  int i;
  float fac1 = 0, fac2 = 0, fac;
  float *ar;

  ar = (float *)malloc(n*sizeof(float));
  for(i=0; i<n1; i++) fac1 += x[i]*x[i];
  for(i=0; i<n2; i++) fac2 += y[i]*y[i];
  fac = sqrt(fac1*fac2);
  if (fac > 1e-12)
  for(i=0; i<n; i++)
    ar[i] = z[i]/fac;
  return ar;
}

float *crscor_0(       float *ar1,
		       float *ar2,
		 const int   n1,
		 const int   n2 ){
   int i, nlen, nwin, wlen, nfft;
   float *ar,*tmp;
   char error[256];
   nlen = n1;
   nwin = 1;
   wlen = nlen;
   nfft = 0;
   tmp = (float *)malloc(4*nlen*sizeof(float));
   crscor(ar1, ar2, nlen,
          nwin, wlen, SAC_RECTANGLE,
          tmp, &nfft, error, 256);
   ar = (float *)malloc((n1+n2-1)*sizeof(float));
   for(i=0; i<=n1-2; i++)
     ar[i] = tmp[nfft-n1+i+1];
   for(i=0; i<=n2-1; i++)
     ar[n1+i-1] = tmp[i];
   free(tmp);
   return ar;
}

