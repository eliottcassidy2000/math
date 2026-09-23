// partial choice E_S (excursion only at even x = r mod 2^J for r in S), print counts every level (q=3,b=1)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int MM=atoi(argv[1]); int J=atoi(argv[2]); char allow[4096]={0}; for(int i=3;i<argc;i++) allow[atoi(argv[i])&((1<<J)-1)]=1;
  const double L2=log(2.0), L3=log(3.0); float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1; float *A=malloc(sizeof(float)*N),*W=malloc(sizeof(float)*N);
    for(uint64_t c=0;c<N;c+=2){ double cont=0; if(m>1){double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w;} A[c]=(float)(-L2+cont); }
    int jmax=(int)((m*L2)/(2*L3))+2;
    for(uint64_t c=0;c<N;c+=2){ double best=A[c]; uint64_t cc=c; double acc=0; for(int j=1;j<=jmax;j++){ if(m<J||!allow[cc&((1<<J)-1)]) break; cc=(9*cc+4)&mask; acc+=2*L3; double v=acc+A[cc]; if(v<best) best=v; } W[c]=(float)best; }
    for(uint64_t c=1;c<N;c+=2){ double w=W[(3*c+1)&mask]; W[c]=(float)(L3+(w<0?w:0.0)); }
    uint64_t nbad=0; for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) nbad++;
    if(m>=14) { printf("S(J=%d) m=%2d exceptional=%llu log2=%.3f\n",J,m,(unsigned long long)nbad,log2((double)nbad)); fflush(stdout); }
    free(A); if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
