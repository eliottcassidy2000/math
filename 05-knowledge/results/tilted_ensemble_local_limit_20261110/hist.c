// hist.c - histogram of H over the tournaments read from stdin (gentourng format). Output "H count".
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
#define MAXN 11
static u64 f[1<<MAXN][MAXN]; static int inm[MAXN];
int main(int argc,char**argv){ int n=atoi(argv[1]); int full=(1<<n)-1; char line[256];
  int HM=200000; u64 *cnt=calloc(HM,sizeof(u64));
  while(fgets(line,sizeof line,stdin)){ if(line[0]!='0'&&line[0]!='1')continue; int k=0;
    for(int i=0;i<n;i++)inm[i]=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(line[k++]=='1')inm[j]|=1<<i; else inm[i]|=1<<j; }
    for(int S=1;S<=full;S++)for(int v=0;v<n;v++){ if(!(S>>v&1))continue; int R=S&~(1<<v);
      if(R==0){f[S][v]=1;continue;} u64 a=0;int m=R&inm[v];while(m){int w=__builtin_ctz(m);m&=m-1;a+=f[R][w];} f[S][v]=a; }
    u64 H=0;for(int v=0;v<n;v++)H+=f[full][v]; if(H<(u64)HM)cnt[H]++; }
  for(int h=0;h<HM;h++)if(cnt[h])printf("%d %llu\n",h,cnt[h]); return 0; }
