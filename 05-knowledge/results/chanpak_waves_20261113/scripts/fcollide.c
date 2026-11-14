// fcollide.c (wave 22): list all collision pairs F(u)=F(u'), 1<=u<u'<=U, of F(u)=908u+9g(u). Usage: ./fcollide U > pairs.txt
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef struct { uint64_t f; uint32_t u; } P;
static int cmp(const void*a,const void*b){const P*x=a,*y=b; if(x->f<y->f)return -1; if(x->f>y->f)return 1; return (x->u<y->u)?-1:(x->u>y->u);}
int main(int argc,char**argv){
  uint32_t U=(uint32_t)strtoull(argv[1],0,10);
  uint64_t *g=malloc((size_t)(U+1)*8); P*p=malloc((size_t)U*sizeof(P));
  g[0]=1;
  for(uint64_t u=1;u<=U;u++){uint64_t v=u>>1; g[u]=(u&1)?g[v]:g[v]+v+(v&1);}
  for(uint64_t u=1;u<=U;u++){p[u-1].f=908*u+9*g[u]; p[u-1].u=(uint32_t)u;}
  qsort(p,U,sizeof(P),cmp);
  uint64_t pairs=0,grp=0; int maxm=1;
  for(uint64_t i=0;i<U;){uint64_t k=i; while(k+1<U&&p[k+1].f==p[i].f)k++;
    int mult=(int)(k-i+1); if(mult>maxm)maxm=mult;
    if(mult>1){grp++; for(uint64_t a=i;a<=k;a++)for(uint64_t b=a+1;b<=k;b++){printf("%u %u\n",p[a].u,p[b].u);pairs++;}}
    i=k+1;}
  fprintf(stderr,"U=%u groups=%llu pairs=%llu maxmult=%d\n",U,(unsigned long long)grp,(unsigned long long)pairs,maxm);
  return 0;}
