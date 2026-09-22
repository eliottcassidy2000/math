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
    uint64_t M=P[r+1];
    V=(float*)malloc(sizeof(float)*M);
    int kmax=(int)(r*L3/L2)+1;
    for(uint64_t c=0;c<M;c++){
      if(c%3==0){ V[c]=INFINITY; continue; }
      int k0 = (c%3==1)?0:1; double best=INFINITY; uint64_t pw=(k0==0)?1:2;
      for(int k=k0;k<=kmax;k+=2){
        uint64_t t=(pw*c)%M; uint64_t y=((t+M-1)%M)/3;
        if(y%3!=0){ double cont=0.0; if(r>1){ double vp=Vprev[y%P[r]]; if(vp<0) cont=vp; } double val=k*L2-L3+cont; if(val<best) best=val; }
        pw=(pw*4)%M;
      }
      V[c]=(float)best;
    }
    if(r==R){ for(uint64_t c=0;c<M;c++) if(c%3 && !(V[c]<0)) printf("%llu %.6f\n",(unsigned long long)c,V[c]); }
    if(Vprev) free(Vprev); Vprev=V;
  }
  return 0;
}
