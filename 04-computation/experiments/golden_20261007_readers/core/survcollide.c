/* survcollide.c: do two distinct descent-survivor words (all prefixes 3^{a_j} > 2^j, 1<=j<=s) of equal length s and
   weight a ever have c(w) = c(w') mod 3^a (i.e. the two classes merge uniformly at time s)?  Exhaustive for s <= SMAX. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int SMAX; static uint64_t P3[64];
typedef struct { uint64_t key; uint32_t a; } rec;
static rec *buf; static size_t nb, cap;
static int surv(int t, int a){ u128 p2=((u128)1)<<t, p3=1; for(int i=0;i<a;i++){ p3*=3; if(p3>p2) return 1;} return p3>p2; }
static void dfs(int t, int a, u128 c, int target){
  if(t>=1 && !surv(t,a)) return;
  if(t==target){ if(nb==cap){cap*=2; buf=realloc(buf,cap*sizeof(rec));} buf[nb].key=(uint64_t)(c % P3[a]); buf[nb].a=a; nb++; return; }
  dfs(t+1,a,c,target);                         /* even step: c unchanged */
  dfs(t+1,a+1,3*c+(((u128)1)<<t),target);      /* odd step: c -> 3c + 2^t */
}
static int cmp(const void*x,const void*y){ const rec*p=x,*q=y; if(p->a!=q->a) return p->a<q->a?-1:1; return p->key<q->key?-1:p->key>q->key; }
int main(int argc,char**argv){
  SMAX=atoi(argv[1]); P3[0]=1; for(int i=1;i<41;i++) P3[i]=P3[i-1]*3;
  cap=1<<20; buf=malloc(cap*sizeof(rec));
  for(int s=1;s<=SMAX;s++){
    nb=0; dfs(0,0,0,s);
    qsort(buf,nb,sizeof(rec),cmp); long coll=0;
    for(size_t i=1;i<nb;i++) if(buf[i].a==buf[i-1].a && buf[i].key==buf[i-1].key) coll++;
    printf("s=%2d survivors %10zu  collisions (pairs of survivors merging at time s) %ld\n",s,nb,coll); fflush(stdout);
  }
}
