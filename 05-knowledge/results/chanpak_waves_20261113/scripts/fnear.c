// fnear.c (wave 22): near-collision counts #{u<u'<=U : |F(u)-F(u')| <= L}, F(u)=908u+9g(u). Usage: ./fnear U
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int cmp(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){
  uint64_t U=strtoull(argv[1],0,10); uint64_t *g=malloc((U+1)*8),*f=malloc(U*8);
  g[0]=1; for(uint64_t u=1;u<=U;u++){uint64_t v=u>>1; g[u]=(u&1)?g[v]:g[v]+v+(v&1);}
  for(uint64_t u=1;u<=U;u++) f[u-1]=908*u+9*g[u];
  qsort(f,U,8,cmp);
  int Ls[]={0,2,8,18,98,998,9998}; int nL=7;
  for(int t=0;t<nL;t++){uint64_t L=Ls[t],cnt=0,j=0; for(uint64_t i=0;i<U;i++){ while(f[i]-f[j]>L) j++; cnt+=i-j; }
    printf("U=%llu L=%5llu pairs=%12llu  pairs/U=%.5f  pairs/(U*(L+1))=%.6f\n",(unsigned long long)U,(unsigned long long)L,(unsigned long long)cnt,(double)cnt/U,(double)cnt/U/(L+1));}
  return 0;}
