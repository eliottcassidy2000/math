// Q2 (reverse E-moves with choice) 3-adic descent game.
// State: residue c mod 3^{r+1}, 3 !| c. Move k>=0 with 2^k c = 1 mod 3; y=(2^k c-1)/3 must be a 3-adic unit.
// V_r(c) = min over legal paths of length<=r of running-min of sum (k*log2 - log3).
// Bad_r = {c : V_r(c) >= 0}  (no multiplicative descent within r steps).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
int main(int argc,char**argv){
  int R=atoi(argv[1]);
  const double L2=log(2.0), L3=log(3.0);
  uint64_t P[40]; P[0]=1; for(int i=1;i<40;i++) P[i]=P[i-1]*3;
  float *Vprev=NULL, *V=NULL;
  for(int r=1;r<=R;r++){
    uint64_t M=P[r+1]; // modulus for level r
    V=(float*)malloc(sizeof(float)*M);
    double kmax_d = r*L3/L2; int kmax=(int)kmax_d+1;
    uint64_t nbad=0, nunit=0;
    for(uint64_t c=0;c<M;c++){
      if(c%3==0){ V[c]=INFINITY; continue; }
      nunit++;
      int k0 = (c%3==1)?0:1;
      double best=INFINITY;
      uint64_t pw = (k0==0)?1:2; // 2^k mod M
      for(int k=k0;k<=kmax;k+=2){
        uint64_t t=(pw*c)%M; // 2^k c mod 3^{r+1}
        // y=(t-1)/3 mod 3^r
        uint64_t y=((t+M-1)%M)/3; // exact since t==1 mod 3
        if(y%3!=0){
          double step=k*L2-L3;
          double cont=0.0;
          if(r>1){ double vp=Vprev[y % P[r]]; if(vp<0) cont=vp; }
          double val=step+cont;
          if(val<best) best=val;
        }
        pw=(pw*4)%M;
      }
      V[c]=(float)best;
      if(!(best<0)) nbad++;
    }
    printf("r=%2d mod 3^%d: units=%llu bad=%llu frac=%.6g\n",r,r+1,(unsigned long long)nunit,(unsigned long long)nbad,(double)nbad/nunit);
    // list bad classes if few
    if(nbad<=40){ printf("   bad:"); for(uint64_t c=0;c<M;c++) if(c%3 && !(V[c]<0)) printf(" %llu",(unsigned long long)c); printf("\n"); }
    fflush(stdout);
    if(Vprev) free(Vprev); Vprev=V;
  }
  return 0;
}
