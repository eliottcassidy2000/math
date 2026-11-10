// labeled.c - exact joint distribution of (number of back arcs, H) over all 2^C(n,2) labelled tournaments on [n].
// Arc {i<j} is "forward" if i->j.  Output: lines "b H count".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 8
static u64 f[1<<MAXN][MAXN];
int main(int argc,char**argv){
  int n=atoi(argv[1]); int m=n*(n-1)/2; int full=(1<<n)-1;
  int pi[32],pj[32],k=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){pi[k]=i;pj[k]=j;k++;}
  int HMAXV=2048; u64 *cnt=calloc((size_t)(m+1)*HMAXV,sizeof(u64));
  for(u64 code=0; code<(1ULL<<m); code++){
    int inm[MAXN]={0}; int b=0;
    for(int e=0;e<m;e++){ if(code>>e&1){ inm[pi[e]]|=1<<pj[e]; b++; } else inm[pj[e]]|=1<<pi[e]; }  // bit set = back arc j->i
    u64 H=0;
    for(int S=1;S<=full;S++) for(int v=0;v<n;v++){ if(!(S>>v&1)) continue; int R=S&~(1<<v);
      if(R==0){f[S][v]=1;continue;} u64 a=0; int mm=R&inm[v]; while(mm){int w=__builtin_ctz(mm);mm&=mm-1;a+=f[R][w];} f[S][v]=a; }
    for(int v=0;v<n;v++) H+=f[full][v];
    cnt[(size_t)b*HMAXV+H]++;
  }
  for(int b=0;b<=m;b++)for(int h=0;h<HMAXV;h++) if(cnt[(size_t)b*HMAXV+h]) printf("%d %d %llu\n",b,h,cnt[(size_t)b*HMAXV+h]);
  return 0;
}
