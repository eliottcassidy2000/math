// Relaxation zoo: at odd x choose any b in Bset (x -> 3x+b), at even x halve; exceptional count mod 2^m.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int MM=atoi(argv[1]); int nb=argc-2; long B[16]; for(int i=0;i<nb;i++) B[i]=atol(argv[2+i]);
  const double L2=log(2.0), L3=log(3.0); float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1; float *W=malloc(sizeof(float)*N);
    for(uint64_t c=0;c<N;c+=2){ double cont=0; if(m>1){double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w;} W[c]=(float)(-L2+cont); }
    for(uint64_t c=1;c<N;c+=2){ double best=INFINITY; for(int i=0;i<nb;i++){ uint64_t y=((uint64_t)(3*c)+(uint64_t)(B[i]%(long)N+N))&mask; double w=W[y]; double v=L3+(w<0?w:0.0); if(v<best) best=v; } W[c]=(float)best; }
    if(m==MM){ uint64_t nbad=0; for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) nbad++; printf("B={"); for(int i=0;i<nb;i++) printf("%s%ld",i?",":"",B[i]); printf("} m=%d exceptional=%llu\n",m,(unsigned long long)nbad); }
    if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
