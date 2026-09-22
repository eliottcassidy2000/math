// Q1 (forward E-moves with choice) 2-adic descent game.
// From odd c: forced 3c+1. From even c: halve (cost -log2, loses a bit) or 9c+4 (cost 2log3, via odd 3c+1).
// W_m(c), c mod 2^m: min over paths of running-min cumulative cost (a log3 - b log2). Bad = W>=0.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
int main(int argc,char**argv){
  int MM=atoi(argv[1]); int dump=argc>2?atoi(argv[2]):0;
  const double L2=log(2.0), L3=log(3.0);
  float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1;
    float *A=(float*)malloc(sizeof(float)*N);
    float *W=(float*)malloc(sizeof(float)*N);
    // A(c) for even c: halve then continue
    for(uint64_t c=0;c<N;c+=2){
      double cont=0.0;
      if(m>1){ double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w; }
      A[c]=(float)(-L2+cont);
    }
    int jmax=(int)((m*L2)/(2*L3))+2;
    for(uint64_t c=0;c<N;c+=2){
      double best=A[c]; uint64_t cc=c;
      for(int j=1;j<=jmax;j++){ cc=(9*cc+4)&mask; double v=2*j*L3+A[cc]; if(v<best) best=v; }
      W[c]=(float)best;
    }
    uint64_t nbad=0;
    for(uint64_t c=1;c<N;c+=2){
      double w=W[(3*c+1)&mask]; double v=L3+(w<0?w:0.0); W[c]=(float)v;
    }
    for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) nbad++;
    printf("m=%2d mod 2^%d: bad=%llu frac=%.3g",m,m,(unsigned long long)nbad,(double)nbad/N);
    if(nbad<=24){ printf("  bad:"); for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) printf(" %lld",(long long)(c>N/2? (int64_t)c-(int64_t)N:(int64_t)c)); }
    printf("\n"); fflush(stdout);
    if(dump && m==MM){ FILE*f=fopen("q1bad_dump.txt","w"); for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) fprintf(f,"%llu %.5f\n",(unsigned long long)c,W[c]); fclose(f);}    
    free(A); if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
