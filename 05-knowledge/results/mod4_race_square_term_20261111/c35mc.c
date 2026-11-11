// c35mc.c - plain Monte Carlo for the parities of the numbers of 3-cycles and 5-cycles of a uniformly random labelled
// tournament on n vertices (n <= 20).  c_3 from the scores; c_5 by enumeration of the 5-cycles through their least
// vertex.  Usage: c35mc n N seed     Output: n N  sum(-1)^c3  sum(-1)^c5  sum(-1)^(c3+c5)
#include <stdio.h>
#include <stdlib.h>
typedef unsigned long long u64; typedef long long i64;
static u64 sm_state;
static inline u64 splitmix(void){ u64 z=(sm_state+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); }
int main(int argc,char**argv){ int n=atoi(argv[1]); u64 N=strtoull(argv[2],0,10); sm_state=strtoull(argv[3],0,10)*0x2545F4914F6CDD1DULL+777;
  unsigned outm[20]; i64 s3=0,s5=0,s35=0;
  for(u64 it=0;it<N;it++){
    for(int i=0;i<n;i++) outm[i]=0; u64 bits=0; int nb=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(nb==0){bits=splitmix();nb=64;} if(bits&1) outm[i]|=1u<<j; else outm[j]|=1u<<i; bits>>=1; nb--; }
    i64 c3=(i64)n*(n-1)*(n-2)/6; for(int i=0;i<n;i++){ int s=__builtin_popcount(outm[i]); c3-=(i64)s*(s-1)/2; }
    u64 c5=0;
    for(int v=0;v<n;v++){ unsigned hi=~((1u<<(v+1))-1) & ((1u<<n)-1);   // vertices > v
      unsigned A=outm[v]&hi;
      for(unsigned a=A;a;a&=a-1){ int x=__builtin_ctz(a); unsigned B=outm[x]&hi&~(1u<<x);
        for(unsigned b=B;b;b&=b-1){ int y=__builtin_ctz(b); unsigned C=outm[y]&hi&~((1u<<x)|(1u<<y));
          for(unsigned c=C;c;c&=c-1){ int z=__builtin_ctz(c); unsigned Dd=outm[z]&hi&~((1u<<x)|(1u<<y)|(1u<<z));
            // last vertex w must beat v
            for(unsigned d=Dd;d;d&=d-1){ int w=__builtin_ctz(d); if(outm[w]>>v&1) c5++; } } } } }
    int p3=(int)(c3&1), p5=(int)(c5&1);
    s3+= p3?-1:1; s5+= p5?-1:1; s35+= (p3^p5)?-1:1;
  }
  printf("%d %llu %lld %lld %lld\n", n, N, s3, s5, s35); return 0; }
