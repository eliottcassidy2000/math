// Generalized forward exceptional-set counter for maps x -> (q x + b) (x odd) / halve (x even),
// with optional excursion at even x: x -> q x + b -> q^2 x + (q+1) b (cost 2 log q), allowed iff x mod 2^J in S.
// usage: excgen q b MM mode [J r1 r2 ...]   mode: 0=no choice, 1=all evens, 2=subset S
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  long q=atol(argv[1]), b=atol(argv[2]); int MM=atoi(argv[3]); int mode=atoi(argv[4]);
  int J=1; char allow[4096]={0};
  if(mode==2){ J=atoi(argv[5]); for(int i=6;i<argc;i++) allow[atoi(argv[i])&((1<<J)-1)]=1; }
  const double L2=log(2.0), LQ=log((double)q);
  float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1;
    float *A=malloc(sizeof(float)*N), *W=malloc(sizeof(float)*N);
    for(uint64_t c=0;c<N;c+=2){ double cont=0; if(m>1){ double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w; } A[c]=(float)(-L2+cont); }
    int jmax=(int)((m*L2)/(2*LQ))+2;
    uint64_t qq=(uint64_t)(q*q)&mask, add=(uint64_t)(((q+1)*b)%(long)N+N)&mask;
    for(uint64_t c=0;c<N;c+=2){
      double best=A[c]; uint64_t cc=c; double acc=0;
      if(mode!=0) for(int j=1;j<=jmax;j++){
        if(mode==2 && (m<J || !allow[cc&((1<<J)-1)])) break;
        cc=(qq*cc+add)&mask; acc+=2*LQ; double v=acc+A[cc]; if(v<best) best=v; }
      W[c]=(float)best;
    }
    uint64_t bb=(uint64_t)((b%(long)N)+N)&mask;
    for(uint64_t c=1;c<N;c+=2){ double w=W[((uint64_t)q*c+bb)&mask]; W[c]=(float)(LQ+(w<0?w:0.0)); }
    uint64_t nbad=0; for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) nbad++;
    if(m==MM || m==MM-4 || m==MM-8) printf("q=%ld b=%ld mode=%d m=%2d exceptional=%llu frac=%.4g log2/m=%.3f\n",q,b,mode,m,(unsigned long long)nbad,(double)nbad/N,nbad?log2((double)nbad)/m:0);
    free(A); if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
