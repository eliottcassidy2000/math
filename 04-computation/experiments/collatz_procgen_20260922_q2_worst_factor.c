// For the backward E-game: V_r(c) = best (most negative) running-min log-multiplier within r steps.
// Report, over classes c mod 3^(r+1) OUTSIDE the hostile neighbourhoods {1 mod 27, 14 mod 27},
// the worst (max) V_r, i.e. the weakest guaranteed descent factor exp(V_r).  Also the histogram.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int R=atoi(argv[1]);
  const double L2=log(2.0), L3=log(3.0);
  uint64_t P[40]; P[0]=1; for(int i=1;i<40;i++) P[i]=P[i-1]*3;
  float *Vprev=NULL,*V=NULL;
  for(int r=1;r<=R;r++){
    uint64_t M=P[r+1]; V=malloc(sizeof(float)*M);
    int kmax=(int)(r*L3/L2)+2;
    for(uint64_t c=0;c<M;c++){
      if(c%3==0){V[c]=INFINITY;continue;}
      int k0=(c%3==1)?0:1; double best=INFINITY; uint64_t pw=(k0==0)?1:2;
      for(int k=k0;k<=kmax;k+=2){ uint64_t t=(pw*c)%M; uint64_t y=((t+M-1)%M)/3;
        if(y%3!=0){ double cont=0; if(r>1){double vp=Vprev[y%P[r]]; if(vp<0) cont=vp;} double val=k*L2-L3+cont; if(val<best) best=val; }
        pw=(pw*4)%M; }
      V[c]=(float)best;
    }
    if(r>=6 && (r%2==0 || r==R)){
      double worst=-1e9; uint64_t wc=0; long hist[8]={0};
      for(uint64_t c=0;c<M;c++){ if(c%3==0) continue; uint64_t m27=c%27; if(m27==1||m27==14) continue;
        double v=V[c]; if(v>worst){worst=v;wc=c;}
        double f=exp(v); int b= f<0.1?0: f<0.2?1: f<0.3?2: f<0.4?3: f<0.5?4: f<0.7?5: f<1.0?6:7; hist[b]++; }
      printf("r=%2d: worst guaranteed factor outside {1,14 mod 27}: exp(V)=%.4f at class %llu (mod 3^%d);  hist[<.1,<.2,<.3,<.4,<.5,<.7,<1,>=1]=",r,exp(worst),(unsigned long long)wc,r+1);
      for(int i=0;i<8;i++) printf("%ld ",hist[i]); printf("\n"); fflush(stdout);
    }
    if(Vprev) free(Vprev); Vprev=V;
  }
  return 0;
}
