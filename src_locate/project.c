/*
 *    code for projecting microseismic sources
 *
 *    Author : Xiaohan Chen
 *
 *    Revision History:
 *      2014.07.07 : initial coding
 *      2014.08.06 : big change 
 *      2014.08.07 : add a direct projection option
 *      2014.09.26 : update to 2.0 version
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

void read_predata();
double distaz(double,double,double,double);
#define PI 3.1415926

int main( int argc, char **argv ){
  int i, x1, x2, y1, y2;
  char label[256];
  float velocity = 3.3;
  
  // read the segment
  for( i=1; i<argc; i++ )
    if ( argv[i][0]=='-' )
      switch ( argv[i][1] ){
        case 'L':
          sscanf( &argv[i][2], "%s", label );
        break;
        case 'R':
          sscanf( &argv[i][2], "%d/%d/%d/%d", &x1, &x2, &y1, &y2 );
        break;
        case 'V':
          sscanf( &argv[i][2], "%f", &velocity );
        break;
        }

  // read previous data
  int num;
  float *tt_c, *dist, *stlo1, *stla1, *stlo2, *stla2, *peak, *coef, *SNR;
  int *SAC_file_1, *SAC_file_2;
  read_predata( label, &num, &tt_c, &dist, &stlo1, &stla1, &stlo2, &stla2, &SAC_file_1, &SAC_file_2, &peak, &coef, &SNR );

  // project
  int k, m;
  float *proj;
  proj = (float *)malloc(sizeof(float)*(x2-x1)*(y2-y1)*100);
  memset(proj,0,sizeof(float)*(x2-x1)*(y2-y1)*100);
  for( i=0; i<num; i++ )
    for(k=1; k<=(x2-x1)*10; k++){
    #pragma omp parallel for firstprivate( k, x1, y1, y2, stla1, stlo1, stla2, stlo2, tt_c, i ) shared (proj)
    for(m=1; m<=(y2-y1)*10; m++){
      float x0 = x1 + 0.1*(k-1) + 0.05;
      float y0 = y1 + 0.1*(m-1) + 0.05;
      float dis0 = (float)distaz( (double) stla1[i],
                                  (double) stlo1[i],
                                  (double) y0,
                                  (double) x0 );
      float dis1 = (float)distaz( (double) stla2[i],
                                  (double) stlo2[i],
                                  (double) y0,
                                  (double) x0 );
      float residual = (dis1-dis0) / velocity - tt_c[i];
      proj[(k-1)*10*(y2-y1)+m-1] += residual * residual;
      }
    }
   
  float ans_lo, ans_la, min_RMS = 999999999, max_RMS;
  for(k=1; k<=(x2-x1)*10; k++)
  for(m=1; m<=(y2-y1)*10; m++){
    float x0 = x1 + 0.1*(k-1) + 0.05;
    float y0 = y1 + 0.1*(m-1) + 0.05;
    float RMS_value = (float) sqrt( proj[(k-1)*10*(y2-y1)+m-1]/(float)num ); 
    if ( RMS_value < min_RMS ){
      ans_lo = x0;
      ans_la = y0;
      min_RMS = RMS_value;
      }
    }
  max_RMS = min_RMS * 2;
  
  // output answer
  FILE *fp;
  char fname[256];
  /* output source location */
  sprintf( fname, "loc_%s", label );
  fp = fopen( fname, "w" );
  fprintf( fp, "%f %f\n%f\n%f\n", ans_lo, ans_la, min_RMS, max_RMS );
  fclose( fp );
  /* output RMS */
  sprintf( fname, "RMS_%s", label );
  fp = fopen( fname, "w" );
  for(k=1; k<=(x2-x1)*10; k++)
  for(m=1; m<=(y2-y1)*10; m++){
    float x0 = x1 + 0.1*(k-1) + 0.05;
    float y0 = y1 + 0.1*(m-1) + 0.05;
    float RMS_value = (float) sqrt( proj[(k-1)*10*(y2-y1)+m-1]/(float)num ); 
    if ( RMS_value <= max_RMS ) fprintf( fp, "%f %f %f\n", x0, y0, RMS_value );
    }
  fclose( fp );
  /* output velocity */
  sprintf( fname, "velo_%s", label );
  fp = fopen( fname, "w" );
  for( i=0; i<num; i++ ){
    float dis0 = (float)distaz( (double) stla1[i],
                                (double) stlo1[i],
                                (double) ans_la,
                                (double) ans_lo );
    float dis1 = (float)distaz( (double) stla2[i],
                                (double) stlo2[i],
                                (double) ans_la,
                                (double) ans_lo );
    float dist_diff = fabs( dis0 - dis1 );
    fprintf( fp, "%f %f\n", fabs(tt_c[i]), dist_diff );
    }
  fclose( fp );

  free(tt_c);
  free(dist);
  free(stlo1);
  free(stla1);
  free(stlo2);
  free(stla2);
  free(peak);
  free(coef);
  free(SNR);
  free(SAC_file_1);
  free(SAC_file_2);
  free(proj);
  return 0;
}

void read_predata( char *label, int *p_num,
                   float **p_tt_c, float **p_dist,
                   float **p_stlo1, float **p_stla1, float **p_stlo2, float **p_stla2,
                   int **p_SAC_file_1, int **p_SAC_file_2,
                   float **p_peak, float **p_coef, float **p_SNR ){
  int i, num;
  FILE *fp;
  char dir[256];

  /* read num of pairs */
  sprintf( dir, "num_picked_pairs_%s", label );
  fp = fopen( dir, "r" );
  fscanf( fp, "%d", &num );
  fclose(fp);
  *p_num = num;
  
  /* malloc */
  float *tt_c = (float *)malloc( num*sizeof(float) ); *p_tt_c = tt_c;
  float *dist = (float *)malloc( num*sizeof(float) ); *p_dist = dist;
  float *stlo1= (float *)malloc( num*sizeof(float) ); *p_stlo1 = stlo1;
  float *stla1= (float *)malloc( num*sizeof(float) ); *p_stla1 = stla1;
  float *stlo2= (float *)malloc( num*sizeof(float) ); *p_stlo2 = stlo2;
  float *stla2= (float *)malloc( num*sizeof(float) ); *p_stla2 = stla2;
  float *peak = (float *)malloc( num*sizeof(float) ); *p_peak = peak;
  float *coef = (float *)malloc( num*sizeof(float) ); *p_coef = coef;
  float *SNR  = (float *)malloc( num*sizeof(float) ); *p_SNR  = SNR;
  int *SAC_file_1 = (int *)malloc( num*sizeof(int) ); *p_SAC_file_1 = SAC_file_1;
  int *SAC_file_2 = (int *)malloc( num*sizeof(int) ); *p_SAC_file_2 = SAC_file_2;
  
  /* read previous answer */
  sprintf( dir, "answer_%s", label );
  fp = fopen( dir, "r" );
  for( i=0; i<num; i++ ){
    //read
    fscanf( fp, "%f %f %f %f %f %f %d %d %f %f %f\n", &tt_c[i], &dist[i], &stlo1[i], &stla1[i], &stlo2[i], &stla2[i], &SAC_file_1[i], &SAC_file_2[i], &peak[i], &coef[i], &SNR[i] );
    }
  fclose(fp);

}

