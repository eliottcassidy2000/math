// Maximal gaps between consecutive attained odd values of one-apex readouts R(U) = a*g(U) + c*(1+U) - s*E(U)
// (E(U)=1 for U even incl. 0, else 0), and of plain g (a=1,c=0,s=0), per dyadic window.
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 g(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((u64)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
int main(int argc,char**argv){
  int K=atoi(argv[1]); i64 a=atoll(argv[2]), c=atoll(argv[3]), s=atoll(argv[4]);
  i64 X=1LL<<K; unsigned char *cov=calloc(X/2+2,1);
  i64 Umax = c>0 ? X/c+2 : 2*X+2;
  for(i64 U=0;U<=Umax;U++){ i64 v=a*g(U)+c*(1+U)-s*((U&1)?0:1); if(v>0 && v<=X && (v&1)) cov[v>>1]=1; }
  printf("R(U) = %lld g(U) + %lld (1+U) - %lld E(U):\n",a,c,s);
  for(int k=4;k<K;k++){ i64 lo=1LL<<k, hi=1LL<<(k+1); i64 prev=-1, mg=0, cnt=0;
    for(i64 n=lo|1;n<hi;n+=2) if(cov[n>>1]){ cnt++; if(prev>=0 && n-prev>mg) mg=n-prev; prev=n; }
    printf("  [2^%d,2^%d): covered %.3f of odds, max gap %lld\n",k,k+1,(double)cnt/((hi-lo)/2),mg); }
  return 0; }
