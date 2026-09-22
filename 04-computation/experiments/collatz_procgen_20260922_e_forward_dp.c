// Partial-choice E_S: Collatz T plus excursion (even c -> 3c+1 -> 9c+4) allowed iff c mod 2^J in S.
// usage: q1partial MM J r1 r2 ...   (residues mod 2^J; must be even)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int MM=atoi(argv[1]); int J=atoi(argv[2]); uint64_t JM=(1ULL<<J)-1;
  char allow[1<<12]={0}; for(int i=3;i<argc;i++) allow[atoi(argv[i])&JM]=1;
  const double L2=log(2.0), L3=log(3.0);
  float *Wprev=NULL;
  for(int m=1;m<=MM;m++){
    uint64_t N=1ULL<<m, mask=N-1;
    float *A=(float*)malloc(sizeof(float)*N), *W=(float*)malloc(sizeof(float)*N);
    for(uint64_t c=0;c<N;c+=2){ double cont=0; if(m>1){ double w=Wprev[(c>>1)&(mask>>1)]; if(w<0) cont=w; } A[c]=(float)(-L2+cont); }
    int jmax=(int)((m*L2)/(2*L3))+2;
    for(uint64_t c=0;c<N;c+=2){
      double best=A[c];
      // excursion chain: only allowed from states whose residue (at the time) is in S; residue mod 2^J known iff m>=J
      uint64_t cc=c; double acc=0; 
      for(int j=1;j<=jmax;j++){
        if(m<J) break; if(!allow[cc&JM]) break;
        cc=(9*cc+4)&mask; acc+=2*L3; double v=acc+A[cc]; if(v<best) best=v;
      }
      W[c]=(float)best;
    }
    for(uint64_t c=1;c<N;c+=2){ double w=W[(3*c+1)&mask]; W[c]=(float)(L3+(w<0?w:0.0)); }
    uint64_t nbad=0; for(uint64_t c=0;c<N;c++) if(!(W[c]<0)) nbad++;
    if(m%2==0 || m==MM) printf("m=%2d bad=%llu log2/m=%.3f\n",m,(unsigned long long)nbad,nbad?log2((double)nbad)/m:0.0);
    free(A); if(Wprev) free(Wprev); Wprev=W;
  }
  return 0;
}
