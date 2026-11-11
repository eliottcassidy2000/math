/* pairspec.c - the pair (H(X), Pi_2(X)) of every tournament class X of order m (code = first field of each input line).
 * Pi_2(X) = sum over vertex sets S of H(X[S]) H(X - S) = number of ordered pairs of disjoint possibly empty paths covering X.
 * usage: pairspec m < cls_m.txt      prints "m code H Pi2" */
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64;
#define MAXV 10
static u64 en[1<<MAXV][MAXV], hs[1<<MAXV];
int main(int argc,char**argv){ int m=atoi(argv[1]); char line[1<<16]; int full=(1<<m)-1;
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); unsigned out[MAXV]={0}; int b=0;
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1) out[i]|=1u<<j; else out[j]|=1u<<i; b++; }
    for(int s=0;s<=full;s++) for(int v=0;v<m;v++) en[s][v]=0;
    for(int v=0;v<m;v++) en[1<<v][v]=1;
    for(int s=1;s<=full;s++) for(int v=0;v<m;v++){ u64 c=en[s][v]; if(c){unsigned nx=out[v]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}} }
    hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int v=0;v<m;v++) t+=en[s][v]; hs[s]=t;}
    u64 p=0; for(int s=0;s<=full;s++) p+=hs[s]*hs[full^s];
    printf("%d %llu %llu %llu\n",m,code,hs[full],p); }
  return 0; }
