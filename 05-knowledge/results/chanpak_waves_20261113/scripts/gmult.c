// Multiplicity of one-apex values over canonical words (U even, or U=0): r(N) = #{U canonical: g(U)=N}, N <= X.
// Reports, per dyadic window, the number of values hit, the total words, the max multiplicity, and sum r^2.
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
static i64 g(i64 U){ if(U==0) return 1; int k=64-__builtin_clzll((u64)U); i64 s=1;
  for(int t=0;t<k;t++){ int bt=(U>>t)&1; if(!bt) s+=U>>(t+1); if(t+1<k && ((U>>(t+1))&1) && !bt) s+=1; } return s; }
int main(int argc,char**argv){
  int K=atoi(argv[1]); i64 X=1LL<<K; unsigned short *r=calloc(X/2+2,sizeof(unsigned short));
  for(i64 U=0;U<=2*X+2;U+=2){ i64 v=g(U); if(v<=X) r[v>>1]++; }
  for(int k=8;k<K;k++){ i64 lo=1LL<<k,hi=1LL<<(k+1); i64 hit=0,words=0,mx=0; double s2=0;
    for(i64 n=lo|1;n<hi;n+=2){ i64 c=r[n>>1]; if(c){hit++; words+=c; s2+=(double)c*c; if(c>mx) mx=c;} }
    double odd=(hi-lo)/2.0;
    printf("  [2^%d,2^%d): hit %.3f of odds, words/odd %.3f, max mult %lld, CS bound words^2/sum r^2 / odds = %.3f\n",
      k,k+1,hit/odd,words/odd,mx,(double)words*words/s2/odd); }
  return 0; }
