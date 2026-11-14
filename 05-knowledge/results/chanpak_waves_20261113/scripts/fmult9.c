// fmult9.c (wave 22): collision census of F(u)=908u+9g(u) on [1,U], one residue class of u mod 9 at a time. Usage: ./fmult9 U
// (F(u)=F(u') forces 917(u'-u) = 9(...), so 9 | u'-u). g via the suffix formula with 14-bit tables.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int cmp(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){
  uint64_t U=strtoull(argv[1],0,10); const int R=14; const uint64_t M=1u<<R;
  uint64_t *gt=malloc(M*8),*rv=malloc(M*8); gt[0]=1;
  for(uint64_t u=1;u<M;u++){uint64_t v=u>>1; gt[u]=(u&1)?gt[v]:gt[v]+v+(v&1);}
  for(uint64_t w=0;w<M;w++){uint64_t r=0; for(int i=0;i<R;i++) if(w>>i&1) r|=1ull<<(R-1-i); rv[w]=r;}
  uint64_t *f=malloc((U/9+2)*8); uint64_t hist[8]={0},pairs=0;
  for(uint64_t cls=0;cls<9;cls++){
    uint64_t n=0;
    for(uint64_t u=(cls?cls:9);u<=U;u+=9){uint64_t a=u>>R,w=u&(M-1),g;
      if(a==0) g=gt[w]; else g=gt[a]+gt[w]-1+a*(M-1-rv[w])+(a&1)*(1-((w>>(R-1))&1));
      f[n++]=908*u+9*g;}
    qsort(f,n,8,cmp);
    for(uint64_t i=0;i<n;){uint64_t k=i; while(k+1<n&&f[k+1]==f[i])k++; uint64_t m=k-i+1; hist[m<7?m:7]++; pairs+=m*(m-1)/2; i=k+1;}
  }
  printf("U=%llu pairs=%llu pairs/U=%.6f\n",(unsigned long long)U,(unsigned long long)pairs,(double)pairs/U);
  for(int m=1;m<8;m++) if(hist[m]) printf("  multiplicity %d: %llu values\n",m,(unsigned long long)hist[m]);
  return 0;}
