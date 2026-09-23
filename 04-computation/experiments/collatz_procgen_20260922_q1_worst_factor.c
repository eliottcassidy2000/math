// Forward E-game: worst guaranteed descent factor over classes mod 2^m OUTSIDE x = -1 mod 2^J.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int MM=atoi(argv[1]);
  const double L2=log(2.0), L3=log(3.0);
  float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1;
    float *A=malloc(sizeof(float)*N),*W=malloc(sizeof(float)*N);
    for(uint64_t c=0;c<N;c+=2){ double cont=0; if(m>1){double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w;} A[c]=(float)(-L2+cont); }
    int jmax=(int)((m*L2)/(2*L3))+2;
    for(uint64_t c=0;c<N;c+=2){ double best=A[c]; uint64_t cc=c; for(int j=1;j<=jmax;j++){ cc=(9*cc+4)&mask; double v=2*j*L3+A[cc]; if(v<best) best=v; } W[c]=(float)best; }
    for(uint64_t c=1;c<N;c+=2){ double w=W[(3*c+1)&mask]; W[c]=(float)(L3+(w<0?w:0.0)); }
    if(m>=10 && m%2==0){
      for(int J=2;J<=8;J+=1){
        double worst=-1e9; uint64_t wc=0; uint64_t nbad=0; uint64_t JM=(1ULL<<J)-1;
        for(uint64_t c=0;c<N;c++){ if(((c+1)&JM)==0) continue; double v=W[c]; if(v>worst){worst=v;wc=c;} if(!(v<0)) nbad++; }
        printf("m=%2d outside -1 mod 2^%d: worst factor exp(W)=%.4f (class %llu), exceptional outside=%llu\n",m,J,exp(worst),(unsigned long long)wc,(unsigned long long)nbad);
      }
      fflush(stdout);
    }
    free(A); if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
